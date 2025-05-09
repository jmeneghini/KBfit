#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#include "xml_handler.h"
#include <functional>
#include <utility>
#include <vector>

/*---------------------------------------------------------------
  Parameters for adaptive bracket + Newton‑secant polish
----------------------------------------------------------------*/
struct AdaptiveBracketConfig {
  double initial_step_percent = 1e-2;
  double x_tol = 1e-8;
  double zero_tol = 1e-6;

  double min_step_percent = 1e-5;
  double max_step_percent = 5e-2;
  double step_scale_limit = 3.0; // |h_new / h_old| <= this

  double plateau_mod2_threshold = 1e-4; // “flat” if |Z|^2 above this
  int plateau_count_before_jump = 4; // consecutive flats before ×2
};


/*---------------------------------------------------------------
  Abstract root finder
----------------------------------------------------------------*/
class RootFinder {
public:
  virtual ~RootFinder() = default;
  virtual bool findRoots(double a, double b, std::vector<double>& roots) = 0;
};

/*---------------------------------------------------------------
  Concrete adaptive‑bracket implementation
----------------------------------------------------------------*/
class AdaptiveBracketRootFinder : public RootFinder {
public:
  using EvalFn = std::function<std::complex<double>(double)>;
  // returns {mod2 , imag}

  static AdaptiveBracketConfig makeConfigFromXML(XMLHandler& xmlin);

  AdaptiveBracketRootFinder(const AdaptiveBracketConfig& params,
                            EvalFn& evaluate);

  AdaptiveBracketRootFinder(XMLHandler& xmlin, EvalFn& evalZ)
      : AdaptiveBracketRootFinder(makeConfigFromXML(xmlin), evalZ) {}

  bool findRoots(double a, double b, std::vector<double>& roots) override;

  std::size_t evalCount() const { return eval_count_; }

private:
  struct Zval {
    std::complex<double> value;   // Ω(x)
    double               mod2;    // |Ω|²   (cached to avoid norm() calls)
  };

  const AdaptiveBracketConfig params_;
  const EvalFn                evaluate_;
  mutable std::size_t         eval_count_ = 0;

  /* ---- one‑call cache ------------------------------------ */
  class EvalCache {
  public:
    explicit EvalCache(AdaptiveBracketRootFinder& rf) : rf_(rf) {}
    const Zval& operator()(double x);

  private:
    AdaptiveBracketRootFinder&        rf_;
    std::unordered_map<double, Zval>  cache_;
  };

  /* ---- helpers ------------------------------------------- */
  double polishSecant(double a, double b, EvalCache& Z) const;
};

#endif // ROOT_FINDER_H