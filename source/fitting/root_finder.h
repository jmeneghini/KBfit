#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#include "xml_handler.h"
#include <bit>
#include <complex>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

/*---------------------------------------------------------------
  Parameters for adaptive bracket + Newton‑secant polish
----------------------------------------------------------------*/
struct AdaptiveBracketConfig {
  double initial_step_percent = 1e-2;
  double x_tol = 1e-8;
  double zero_tol = 1e-6;

  double min_step_percent = 1e-3;
  double max_step_percent = 1e-2;
  double step_scale_limit = 3.0; // |h_new / h_old| <= this

  double plateau_mod2_threshold = 0.75; // “flat” if |Z|² above this
  int plateau_count_before_jump = 3;    // consecutive flats before ×2
};

/*---------------------------------------------------------------
  Abstract root finder interface
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

  /*---- construction -----------------------------------------------------*/
  static AdaptiveBracketConfig makeConfigFromXML(XMLHandler& xmlin);

  AdaptiveBracketRootFinder(const AdaptiveBracketConfig& params,
                            EvalFn& evaluate);

  AdaptiveBracketRootFinder(XMLHandler& xmlin, EvalFn& evalZ)
      : AdaptiveBracketRootFinder(makeConfigFromXML(xmlin), evalZ) {}

  /*---- API --------------------------------------------------------------*/
  bool findRoots(double a, double b, std::vector<double>& roots) override;

  std::size_t evalCount() const { return eval_count_; }

private:
  /*---- cached Ω evaluation ---------------------------------------------*/
  struct Zval {
    std::complex<double> value; // Ω(x)
    double mod2;                // |Ω|²
  };

  /* bit‑wise hashing so identical IEEE‑754 bits share one cache entry */
  struct KeyHash {
    std::size_t operator()(double x) const noexcept {
      return std::hash<std::uint64_t>{}(std::bit_cast<std::uint64_t>(x));
    }
  };
  struct KeyEq {
    bool operator()(double a, double b) const noexcept {
      return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
    }
  };

  class EvalCache {
  public:
    explicit EvalCache(AdaptiveBracketRootFinder& rf) : rf_(rf) {}
    const Zval& operator()(double x);

  private:
    AdaptiveBracketRootFinder& rf_;
    std::unordered_map<double, Zval, KeyHash, KeyEq> cache_;
  };

  /*---- helpers ---------------------------------------------------------*/
  double polishSecant(double a, double b, EvalCache& Z) const;

  /*---- data ------------------------------------------------------------*/
  const AdaptiveBracketConfig params_;
  const EvalFn evaluate_;
  mutable std::size_t eval_count_ = 0;
};

#endif /* ROOT_FINDER_H */