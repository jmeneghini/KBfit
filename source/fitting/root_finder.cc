#include "root_finder.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

// caches evaluated values of the function, keeping count of # of evaluations
const AdaptiveBracketRootFinder::Zval&
AdaptiveBracketRootFinder::EvalCache::operator()(double x) {
  auto [it, inserted] = cache_.emplace(x, Zval{});
  if (inserted) {
    ++rf_.eval_count_;
    auto val = rf_.evaluate_(x);
    it->second.value = val;
    it->second.mod2 = std::norm(val);
  }
  return it->second;
}

// construct from config from XML
AdaptiveBracketConfig
AdaptiveBracketRootFinder::makeConfigFromXML(XMLHandler& xmlin) {
  AdaptiveBracketConfig p;
  xmlreadifchild(xmlin, "InitialStepPercent", p.initial_step_percent);
  xmlreadifchild(xmlin, "AbsXTolerance", p.x_tol);
  xmlreadifchild(xmlin, "AbsResidualTolerance", p.zero_tol);
  xmlreadifchild(xmlin, "MinStepPercent", p.min_step_percent);
  xmlreadifchild(xmlin, "MaxStepPercent", p.max_step_percent);
  xmlreadifchild(xmlin, "StepScaleLimit", p.step_scale_limit);
  xmlreadifchild(xmlin, "PlateauMod2Threshold", p.plateau_mod2_threshold);
  xmlreadifchild(xmlin, "PlateauCountBeforeJump", p.plateau_count_before_jump);
  return p;
}

// ctor
AdaptiveBracketRootFinder::AdaptiveBracketRootFinder(
    const AdaptiveBracketConfig& p, EvalFn& evaluate)
    : params_(p), evaluate_(evaluate) {}

// polishSecant – secant + midpoint fallback (simplified and faster)
double AdaptiveBracketRootFinder::polishSecant(double a, double b,
                                               EvalCache& Z) const {
  constexpr std::size_t MAX_IT = 20;
  const double xtol = params_.x_tol;
  const double ftol = params_.zero_tol;

  double x0 = a, f0 = Z(a).value.imag();
  if (std::abs(f0) < ftol) return x0;
  
  double x1 = b, f1 = Z(b).value.imag();
  if (std::abs(f1) < ftol) return x1;
  
  // Simple robust check - if no clear sign change, use bisection
  if (f0 * f1 > 0) {
    // Try bisection once to see if there's actually a root
    double xmid = 0.5 * (a + b);
    double fmid = Z(xmid).value.imag();
    if (std::abs(fmid) < ftol) return xmid;
    
    // If still no sign change after bisection, it's truly invalid
    if (f0 * fmid > 0 && f1 * fmid > 0) {
      throw std::runtime_error("invalid bracket");
    }
    // Use the valid sub-bracket
    if (f0 * fmid <= 0) { x1 = xmid; f1 = fmid; }
    else { x0 = xmid; f0 = fmid; }
  }

  // Secant method with bisection fallback
  for (std::size_t it = 0; it < MAX_IT; ++it) {
    double denom = f1 - f0;
    if (std::abs(denom) < 1e-14 * std::max(std::abs(f0), std::abs(f1)))
      break;

    double x2 = x1 - f1 * (x1 - x0) / denom;
    if (x2 <= a || x2 >= b) x2 = 0.5 * (x0 + x1);

    double f2 = Z(x2).value.imag();
    if (std::abs(f2) < ftol || std::abs(x2 - x1) < xtol) return x2;

    // Update bracket maintaining sign change
    if (f1 * f2 < 0) { x0 = x1; f0 = f1; }
    else { x0 = x2; f0 = f2; }
    x1 = x2; f1 = f2;
  }
  return 0.5 * (x0 + x1);
}

// More robust bracket validation
inline bool isValidBracket(double f_left, double f_right, double zero_tol) {
  // Allow for numerical precision: if either value is very small, consider it valid
  const double precision_factor = 1e-12;
  double threshold = zero_tol * precision_factor;
  
  return (f_left * f_right <= 0) || 
         (std::abs(f_left) < threshold) || 
         (std::abs(f_right) < threshold);
}

// findRoots – main driver (optimized)
bool AdaptiveBracketRootFinder::findRoots(double a, double b,
                                          std::vector<double>& roots) {
  roots.clear();
  eval_count_ = 0;
  EvalCache Z(*this);

  const double L = b - a;
  const double h0 = L * params_.initial_step_percent;
  const double h_min = L * params_.min_step_percent;
  const double h_max = L * params_.max_step_percent;
  const double eps = std::numeric_limits<double>::epsilon();

  std::deque<std::pair<double, double>> brackets;

  // Adaptive scan with optimized bracket collection
  double x = a;
  double h = std::min(h0, L);
  
  auto Zv = Z(x);
  double im_prev = Zv.value.imag();
  double mod2_prev = Zv.mod2;
  int flat_cnt = 0;

  while (x < b - eps) {
    double x_next = std::min(x + h, b);
    auto Zn = Z(x_next);
    double im_curr = Zn.value.imag();
    double mod2_curr = Zn.mod2;

    // Primary bracket detection: direct sign flip
    if (im_prev * im_curr < 0) {
      brackets.emplace_back(x, x_next);
    }
    // Secondary: near-zero magnitude with robust validation
    else if (mod2_curr < params_.zero_tol) {
      double left_b = std::max(a, x - h_max);
      double right_b = std::min(b, x_next + h_max);
      
      // Only evaluate if bounds are reasonable and different
      if (right_b - left_b > params_.x_tol) {
        double f_left = Z(left_b).value.imag();
        double f_right = Z(right_b).value.imag();
        
        if (isValidBracket(f_left, f_right, params_.zero_tol)) {
          brackets.emplace_back(left_b, right_b);
        }
      }
    }
    // large imaginary jump detection
    else {
      double im_diff = std::abs(im_curr - im_prev);
      double jump_threshold = 0.3 * std::sqrt(mod2_prev);
      
      if (im_diff > jump_threshold && mod2_prev > params_.zero_tol) {
        double x_mid = 0.5 * (x + x_next);
        auto Zm = Z(x_mid);
        double im_mid = Zm.value.imag();
        
        // Check sub-intervals for sign changes
        if (im_prev * im_mid < 0) brackets.emplace_back(x, x_mid);
        if (im_mid * im_curr < 0) brackets.emplace_back(x_mid, x_next);
      }
    }

    // step size control
    if (mod2_curr > params_.plateau_mod2_threshold) {
      h = h_max;
      flat_cnt = 0;
    } else {
      // Fast ratio calculation avoiding sqrt when possible
      double ratio;
      if (mod2_prev > eps) {
        ratio = mod2_curr / mod2_prev;
        ratio = (ratio > 1.0) ? std::sqrt(ratio) : 1.0 / std::sqrt(1.0 / ratio);
      } else {
        ratio = 1.0;
      }
      
      ratio = std::clamp(ratio, 1.0 / params_.step_scale_limit, 
                         params_.step_scale_limit);
      h = std::clamp(h * ratio, h_min, h_max);
      
      // Plateau jump logic
      if (ratio > 0.9 && ratio < 1.1) {
        if (++flat_cnt >= params_.plateau_count_before_jump) {
          h = std::min(2.0 * h, h_max);
          flat_cnt = 0;
        }
      } else {
        flat_cnt = 0;
      }
    }

    x = x_next;
    im_prev = im_curr;
    mod2_prev = mod2_curr;
  }

  // Polish brackets and handle failures with efficient DFS
  roots.reserve(brackets.size()); // Pre-allocate likely size
  
  for (auto [lft, rgt] : brackets) {
    try {
      double root = polishSecant(lft, rgt, Z);
      if (Z(root).mod2 < params_.zero_tol) {
        roots.push_back(root);
        continue;
      }
    } catch (const std::runtime_error&) {
      // Polish failed, use DFS fallback
    }

    // DFS fallback with iterative approach (faster than recursion)
    std::deque<std::pair<double, double>> stack;
    stack.emplace_back(lft, rgt);
    
    while (!stack.empty()) {
      auto [l, r] = stack.back();
      stack.pop_back();
      
      if (r - l <= params_.x_tol) {
        double xm = 0.5 * (l + r);
        if (Z(xm).mod2 < params_.zero_tol) {
          roots.push_back(xm);
        }
        continue;
      }
      
      double mid = 0.5 * (l + r);
      double fim = Z(mid).value.imag();
      double fl = Z(l).value.imag();
      double fr = Z(r).value.imag();
      
      if (fl * fim <= 0) stack.emplace_back(l, mid);
      if (fim * fr <= 0) stack.emplace_back(mid, r);
    }
  }

  // Efficient deduplication and sorting
  if (roots.size() > 1) {
    std::sort(roots.begin(), roots.end());
    auto new_end = std::unique(roots.begin(), roots.end(),
                               [this](double a, double b) {
                                 return std::abs(a - b) < params_.x_tol;
                               });
    roots.erase(new_end, roots.end());
  }
  
  roots.shrink_to_fit();
  return !roots.empty();
}
