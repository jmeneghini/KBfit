#include "root_finder.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

/*---------------------------------------------------------------
  caches evaluated values of the function, keeping count of # of evaluations
----------------------------------------------------------------*/
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

/*---------------------------------------------------------------
  construct from config from XML
----------------------------------------------------------------*/
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

/*---------------------------------------------------------------
  ctor
----------------------------------------------------------------*/
AdaptiveBracketRootFinder::AdaptiveBracketRootFinder(
    const AdaptiveBracketConfig& p, EvalFn& evaluate)
    : params_(p), evaluate_(evaluate) {}

/*---------------------------------------------------------------
  polishSecant – secant + midpoint fallback
----------------------------------------------------------------*/
double AdaptiveBracketRootFinder::polishSecant(double a, double b,
                                               EvalCache& Z) const {
  constexpr std::size_t MAX_IT = 20;
  const double xtol = params_.x_tol;
  const double ftol = params_.zero_tol;

  double x0 = a, f0 = Z(a).value.imag();
  if (std::abs(f0) < ftol)
    return x0;
  double x1 = b, f1 = Z(b).value.imag();
  if (std::abs(f1) < ftol)
    return x1;
  if (f0 * f1 > 0) {
    throw std::runtime_error("invalid bracket");
  }

  for (std::size_t it = 0; it < MAX_IT; ++it) {
    double denom = f1 - f0;
    if (std::abs(denom) < 1e-14 * std::max(std::abs(f0), std::abs(f1)))
      break;

    double x2 = x1 - f1 * (x1 - x0) / denom;
    if (x2 <= a || x2 >= b)
      x2 = 0.5 * (x0 + x1);

    double f2 = Z(x2).value.imag();
    if (std::abs(f2) < ftol || std::abs(x2 - x1) < xtol)
      return x2;

    (f1 * f2 < 0) ? (x0 = x1, f0 = f1) : (x0 = x2, f0 = f2);
    x1 = x2;
    f1 = f2;
  }
  return 0.5 * (x0 + x1); // fallback midpoint
}

/*---------------------------------------------------------------
  findRoots – main driver
----------------------------------------------------------------*/
bool AdaptiveBracketRootFinder::findRoots(double a, double b,
                                          std::vector<double>& roots) {
  roots.clear();
  eval_count_ = 0;
  EvalCache Z(*this);

  const double L = b - a;
  const double h0 = L * params_.initial_step_percent;
  const double h_min = L * params_.min_step_percent;
  const double h_max = L * params_.max_step_percent;

  std::deque<std::pair<double, double>> brackets;

  /* ---- adaptive scan to collect brackets ----------------------------- */
  double x = a;
  double h = std::min(h0, b - a);
  auto Zv = Z(x);
  double im_prev = Zv.value.imag();
  double mod2_prev = Zv.mod2;
  int flat_cnt = 0;

  while (x < b - std::numeric_limits<double>::epsilon()) {
    double x_next = std::min(x + h, b);
    auto Zn = Z(x_next);

    bool check_for_im_jump = true;
    if (im_prev * Zn.value.imag() < 0) { // sign flip
      brackets.emplace_back(x, x_next);
      check_for_im_jump = false;
    } else if (Zn.mod2 < params_.zero_tol) { // mod2 grazes 0
      brackets.emplace_back(x - h_max,
                            x_next + h_max); // give an extended bracket
      check_for_im_jump = false;
    }

    if (check_for_im_jump) {
      double im_diff = abs(Zn.value.imag() - im_prev);
      double jump_thres = 0.3 * mod2_prev;
      if (im_diff > jump_thres) {
        double x_mid = 0.5 * (x + x_next);
        auto Zm = Z(x_mid);
        if (im_prev * Zm.value.imag() < 0 || mod2_prev < params_.zero_tol ||
            Zm.mod2 < params_.zero_tol)
          brackets.emplace_back(x, x_mid);
        if (Zm.value.imag() * Zn.value.imag() < 0 ||
            Zm.mod2 < params_.zero_tol || Zn.mod2 < params_.zero_tol)
          brackets.emplace_back(x_mid, x_next);
      }
    }

    // step‑size control
    double h_new;
    if (Zn.mod2 > params_.plateau_mod2_threshold) {
      h_new = h_max;
      flat_cnt = 0;
    } else {
      double ratio = std::sqrt(Zn.mod2 / std::max(mod2_prev, 1e-300));
      ratio = std::clamp(ratio, 1.0 / params_.step_scale_limit,
                         params_.step_scale_limit);
      h_new = std::clamp(h * ratio, h_min, h_max);

      if (0.9 < ratio && ratio < 1.1) {
        if (++flat_cnt >= params_.plateau_count_before_jump) {
          h_new = std::min(2 * h_new, h_max);
          flat_cnt = 0;
        }
      } else
        flat_cnt = 0;
    }

    x = x_next;
    im_prev = Zn.value.imag();
    mod2_prev = Zn.mod2;
    h = h_new;
  }

  /* ---- polish brackets & DFS fallback -------------------------------- */
  for (auto [lft, rgt] : brackets) {
    double root = polishSecant(lft, rgt, Z);
    if (Z(root).mod2 < params_.zero_tol) {
      roots.push_back(root);
      continue;
    }

    std::deque<std::pair<double, double>> stk{{lft, rgt}};
    while (!stk.empty()) {
      auto [l, r] = stk.back();
      stk.pop_back();
      if (r - l <= params_.x_tol) {
        double xm = 0.5 * (l + r);
        if (Z(xm).mod2 < params_.zero_tol)
          roots.push_back(xm);
        continue;
      }
      double mid = 0.5 * (l + r);
      double fim = Z(mid).value.imag();
      if (Z(l).value.imag() * fim <= 0)
        stk.emplace_back(l, mid);
      else if (fim * Z(r).value.imag() <= 0)
        stk.emplace_back(mid, r);
    }
  }

  /* ---- deduplicate --------------------------------------------------- */
  std::sort(roots.begin(), roots.end());
  roots.erase(std::unique(roots.begin(), roots.end(),
                          [&](double p, double q) {
                            return std::abs(p - q) < params_.x_tol;
                          }),
              roots.end());
  /* ----- resize to remove unused capacity ------------------------- */
  roots.shrink_to_fit();
  return !roots.empty();
}