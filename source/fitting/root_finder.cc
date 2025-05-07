/*======================  adaptive_bracket.cpp  =====================*/
#include "root_finder.h"
#include <deque>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <stdexcept>

/*---------------------------------------------------------------
    EvalCache implementation
----------------------------------------------------------------*/
const AdaptiveBracketRootFinder::Zval&
AdaptiveBracketRootFinder::EvalCache::operator()(double x)
{
    auto [it, inserted] = cache_.emplace(x, Zval{});
    if (inserted) {
        ++rf_.n_eval_unique_;
        auto pr          = rf_.evalZ_(x);   // {mod2, imag}
        it->second.mod2  = pr.first;
        it->second.im    = pr.second;
    }
    return it->second;
}

/*---------------------------------------------------------------
   ctor
----------------------------------------------------------------*/
AdaptiveBracketRootFinder::AdaptiveBracketRootFinder(
        const AdaptiveBracketParameters& p,
        EvalFn evalZ)
    : params_(p), evalZ_(std::move(evalZ))
{}

/*---------------------------------------------------------------
   polishSecant – safeguarded Newton/secant in [a,b]
----------------------------------------------------------------*/
double AdaptiveBracketRootFinder::polishSecant(double a, double b,
                                               EvalCache&   Z) const
{
    const std::size_t MAX_IT = 20;
    const double xtol = params_.eps_fine;
    const double ftol = params_.tol_zero;

    double fa = Z(a).im;
    if (std::abs(fa) < ftol) return a;
    double fb = Z(b).im;
    if (std::abs(fb) < ftol) return b;
    if (fa * fb > 0) throw std::runtime_error("invalid bracket");

    double x0 = a, f0 = fa;
    double x1 = b, f1 = fb;

    for (std::size_t it = 0; it < MAX_IT; ++it) {
        double denom = f1 - f0;
        if (denom == 0.0) break;
        double x2 = x1 - f1 * (x1 - x0) / denom;
        if (x2 <= a || x2 >= b) x2 = 0.5 * (a + b);

        double f2 = Z(x2).im;
        if (std::abs(f2) < ftol || std::abs(x2 - x1) < xtol) return x2;

        if (f1 * f2 < 0) { a = x1; fa = f1; b = x2; fb = f2; }
        else             { a = x2; fa = f2; }

        x0 = a; f0 = fa; x1 = b; f1 = fb;
    }
    return 0.5 * (a + b);              // bisection fallback
}

/*---------------------------------------------------------------
   main algorithm
----------------------------------------------------------------*/
bool AdaptiveBracketRootFinder::findRoots(double a, double b,
                                          std::vector<double>& roots)
{
    roots.clear();
    n_eval_unique_ = 0;
    EvalCache Z(*this);

    /* ---------- adaptive coarse scan ------------------------ */
    std::deque<std::pair<double,double>> brackets;

    double x = a;
    double h = std::min(params_.starting_step, b - a);
    auto   Zv = Z(x);
    double im_prev = Zv.im;
    double mod2_prev = Zv.mod2;
    int    flat_cnt = 0;

    while (x < b - std::numeric_limits<double>::epsilon())
    {
        double x_next = std::min(x + h, b);
        auto   Zn     = Z(x_next);

        if (im_prev * Zn.im < 0 || mod2_prev < params_.tol_zero
                                || Zn.mod2     < params_.tol_zero)
            brackets.emplace_back(x, x_next);

        /* ---- adaptive step update -------------------------- */
        double h_new;
        if (Zn.mod2 > params_.Z_flat) {
            h_new = params_.h_max; flat_cnt = 0;
        } else {
            double ratio = std::sqrt(Zn.mod2 / mod2_prev);
            ratio = std::clamp(ratio, 1.0/params_.beta, params_.beta);
            h_new = std::clamp(h * ratio, params_.h_min, params_.h_max);

            if (0.9 < ratio && ratio < 1.1) {
                if (++flat_cnt >= params_.m_flat) {
                    h_new  = std::min(2*h_new, params_.h_max);
                    flat_cnt = 0;
                }
            } else flat_cnt = 0;
        }

        x = x_next; im_prev = Zn.im; mod2_prev = Zn.mod2; h = h_new;
    }

    /* ---------- fine polish --------------------------------- */
    for (auto [lft, rgt] : brackets)
    {
        double root = polishSecant(lft, rgt, Z);
        if (Z(root).mod2 < params_.tol_zero) { roots.push_back(root); continue; }


        std::deque<std::pair<double,double>> stack{{lft, rgt}};
        while (!stack.empty()) {
            auto [l, r] = stack.back(); stack.pop_back();
            if (r - l <= params_.eps_fine) {
                double xm = 0.5*(l+r);
                if (Z(xm).mod2 < params_.tol_zero) roots.push_back(xm);
                continue;
            }
            double mid = 0.5*(l+r);
            if (Z(l).im  * Z(mid).im <= 0) stack.emplace_back(l, mid);
            if (Z(mid).im* Z(r).im  <= 0) stack.emplace_back(mid, r);
        }
    }

    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end(),
                            [&](double p,double q){
                                return std::abs(p-q) < params_.eps_fine;
                            }),
                roots.end());
    return true;
}

/*----------------------------  EOF  ------------------------------*/