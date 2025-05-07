/*=========================  root_finder.h  =========================*/
#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#include <vector>
#include <functional>
#include <utility>   // std::pair

/*---------------------------------------------------------------
  Abstract “parameter pack”
----------------------------------------------------------------*/
class RootFinderParameters {
public:
    virtual ~RootFinderParameters() = default;
};

/*---------------------------------------------------------------
  Parameters for adaptive bracket + Newton‑secant polish
----------------------------------------------------------------*/
class AdaptiveBracketParameters : public RootFinderParameters {
public:
    double starting_step = 1e-2;
    double eps_fine      = 1e-6;
    double tol_zero      = 1e-10;

    double h_min         = 1e-4;
    double h_max         = 5e-2;
    double beta          = 3.0;
    double Z_flat        = 1e-4;
    int    m_flat        = 4;
};

/*---------------------------------------------------------------
  Abstract root finder
----------------------------------------------------------------*/
class RootFinder {
public:
    virtual ~RootFinder() = default;
    virtual bool findRoots(double a,
                           double b,
                           std::vector<double>& roots) = 0;
};

/*---------------------------------------------------------------
  Concrete adaptive‑bracket implementation
----------------------------------------------------------------*/
class AdaptiveBracketRootFinder : public RootFinder {
public:
    using EvalFn = std::function<std::pair<double,double>(double)>;
                    // returns {mod2 , imag}

    AdaptiveBracketRootFinder(const AdaptiveBracketParameters& params,
                              EvalFn evalZ);

    bool findRoots(double a,
                   double b,
                   std::vector<double>& roots) override;

    std::size_t evalCount() const { return n_eval_unique_; }

private:
    struct Zval { double mod2; double im; };

    const AdaptiveBracketParameters params_;
    const EvalFn evalZ_;
    mutable std::size_t n_eval_unique_ = 0;

    /* ---- one‑call cache ------------------------------------ */
    class EvalCache {
    public:
        explicit EvalCache(AdaptiveBracketRootFinder& rf) : rf_(rf) {}
        const Zval& operator()(double x);
    private:
        AdaptiveBracketRootFinder& rf_;
        std::unordered_map<double, Zval> cache_;
    };

    /* ---- helpers ------------------------------------------- */
    double polishSecant(double a, double b, EvalCache& Z) const;
};

#endif // ROOT_FINDER_H