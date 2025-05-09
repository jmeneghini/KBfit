// adaptive_bracket.cpp
#include "root_finder.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>
#include <unordered_map>

const AdaptiveBracketRootFinder::Zval&
AdaptiveBracketRootFinder::EvalCache::operator()(double x) {
  auto [it, inserted] = cache_.emplace(x, Zval{});
  if (inserted) {
    ++rf_.eval_count_;
    std::complex<double> val = rf_.evaluate_(x);
    it->second.value = val;
    it->second.mod2  = std::norm(val);
  }
  return it->second;
}
/*---------------------------------------------------------------
   makeAdaptiveBracketConfigFromXML
----------------------------------------------------------------*/
AdaptiveBracketConfig
AdaptiveBracketRootFinder::makeConfigFromXML(XMLHandler& xmlin) {
  AdaptiveBracketConfig p;
  xmlreadifchild(xmlin, "InitialStepSize", p.initial_step_size);
  xmlreadifchild(xmlin, "AbsXTolerance", p.x_tol);
  xmlreadifchild(xmlin, "AbsResidualTolerance", p.zero_tol);
  xmlreadifchild(xmlin, "MinStepSize", p.min_step_size);
  xmlreadifchild(xmlin, "MaxStepSize", p.max_step_size);
  xmlreadifchild(xmlin, "StepScaleLimit", p.step_scale_limit);
  xmlreadifchild(xmlin, "PlateauMod2Threshold", p.plateau_mod2_threshold);
  xmlreadifchild(xmlin, "PlateauCountBeforeJump", p.plateau_count_before_jump);
  return p;
}

AdaptiveBracketRootFinder::AdaptiveBracketRootFinder(
    const AdaptiveBracketConfig& p, EvalFn& evaluate)
    : params_(p), evaluate_(evaluate) {}

double AdaptiveBracketRootFinder::polishSecant(double a, double b,
                                               EvalCache& Z) const {
  const std::size_t MAX_IT = 20;
  const double xtol = params_.x_tol;
  const double ftol = params_.zero_tol;

  double fa = Z(a).value.imag();
  if (std::abs(fa) < ftol) return a;
  double fb = Z(b).value.imag();
  if (std::abs(fb) < ftol) return b;
  if (fa * fb > 0) throw std::runtime_error("invalid bracket");

  double x0 = a, f0 = fa;
  double x1 = b, f1 = fb;

  for (std::size_t it = 0; it < MAX_IT; ++it) {
    double denom = f1 - f0;
    if (denom == 0.0) break;
    double x2 = x1 - f1 * (x1 - x0) / denom;
    if (x2 <= a || x2 >= b) x2 = 0.5 * (a + b);

    double f2 = Z(x2).value.imag();
    if (std::abs(f2) < ftol || std::abs(x2 - x1) < xtol) return x2;

    if (f1 * f2 < 0) {
      a = x1; fa = f1; b = x2; fb = f2;
    } else {
      a = x2; fa = f2;
    }

    x0 = a; f0 = fa; x1 = b; f1 = fb;
  }
  return 0.5 * (a + b);
}

bool AdaptiveBracketRootFinder::findRoots(double a, double b,
                                          std::vector<double>& roots) {
  roots.clear();
  eval_count_ = 0;
  EvalCache Z(*this);

  std::deque<std::pair<double, double>> brackets;

  double x = a;
  double h = std::min(params_.initial_step_size, b - a);
  auto Zv = Z(x);
  double im_prev  = Zv.value.imag();
  double mod2_prev= Zv.mod2;
  int flat_cnt    = 0;

  while (x < b - std::numeric_limits<double>::epsilon()) {
    double x_next = std::min(x + h, b);
    auto Zn = Z(x_next);

    if (im_prev * Zn.value.imag() < 0 ||
        mod2_prev < params_.zero_tol ||
        Zn.mod2   < params_.zero_tol)
      brackets.emplace_back(x, x_next);

    double h_new;
    if (Zn.mod2 > params_.plateau_mod2_threshold) {
      h_new = params_.max_step_size;
      flat_cnt = 0;
    } else {
      double ratio = std::sqrt(Zn.mod2 / mod2_prev);
      ratio = std::clamp(ratio,
                         1.0 / params_.step_scale_limit,
                         params_.step_scale_limit);
      h_new = std::clamp(h * ratio,
                         params_.min_step_size,
                         params_.max_step_size);

      if (0.9 < ratio && ratio < 1.1) {
        if (++flat_cnt >= params_.plateau_count_before_jump) {
          h_new = std::min(2 * h_new, params_.max_step_size);
          flat_cnt = 0;
        }
      } else {
        flat_cnt = 0;
      }
    }

    x = x_next;
    im_prev  = Zn.value.imag();
    mod2_prev= Zn.mod2;
    h        = h_new;
  }

  for (auto [lft, rgt] : brackets) {
    double root = polishSecant(lft, rgt, Z);
    if (Z(root).mod2 < params_.zero_tol) {
      roots.push_back(root);
      continue;
    }

    std::deque<std::pair<double, double>> stack{{lft, rgt}};
    while (!stack.empty()) {
      auto [l, r] = stack.back();
      stack.pop_back();
      if (r - l <= params_.x_tol) {
        double xm = 0.5 * (l + r);
        if (Z(xm).mod2 < params_.zero_tol) roots.push_back(xm);
        continue;
      }
      double mid = 0.5 * (l + r);
      if (Z(l).value.imag()  * Z(mid).value.imag() <= 0) stack.emplace_back(l, mid);
      if (Z(mid).value.imag()* Z(r).value.imag() <= 0) stack.emplace_back(mid, r);
    }
  }

  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end(),
                          [&](double p, double q) {
                            return std::abs(p - q) < params_.x_tol;
                          }),
              roots.end());
  return true;
}