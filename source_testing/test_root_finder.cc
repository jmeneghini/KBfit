
#ifndef M_PI // Define M_PI if not defined (e.g., on MSVC without _USE_MATH_DEFINES)
    #define M_PI 3.14159265358979323846
#endif

#include "doctest.h"
#include "root_finder.h" // Your header file
#include "helper.h"
#include <cmath>        // For M_PI, sin, cos, pow, sqrt, abs, round
#include <vector>
#include <string>
#include <sstream>      // For ostringstream
#include <algorithm>    // For std::sort
#include <iomanip>      // For std::setprecision
#include <functional>   // For std::function
#include <complex>      // For std::complex

using namespace TestHelper;

// Helper to convert vector to string for Doctest INFO messages
template <typename T>
std::string vec_to_string(const std::vector<T>& vec) {
    std::ostringstream oss;
    oss << "[";
    bool first = true;
    for (const auto& val : vec) {
        if (!first) {
            oss << ", ";
        }
        oss << std::fixed << std::setprecision(10) << val; // Formatting within ostringstream is fine
        first = false;
    }
    oss << "]";
    return oss.str();
}

// Helper function to run a root-finding test and check results
void check_roots(const std::string& test_name,
                 AdaptiveBracketRootFinder::EvalFn func_to_test,
                 double a, double b,
                 const std::vector<double>& expected_roots,
                 const AdaptiveBracketConfig& config = AdaptiveBracketConfig(), // Default config
                 double root_comparison_tol = 1e-7, // Tolerance for comparing root values
                 bool check_eval_count = false, // Optional: check eval count (for specific performance checks)
                 size_t max_expected_evals = 0) { // Renamed from check_eval_count to avoid conflict

    AdaptiveBracketRootFinder::EvalFn local_func = func_to_test; // EvalFn taken by non-const ref by finder
    AdaptiveBracketRootFinder finder(config, local_func);
    std::vector<double> found_roots;

    // Pre-format complex strings for INFO macro
    std::ostringstream oss_interval_info;
    oss_interval_info << "Interval: [" << std::fixed << std::setprecision(10) << a
                      << ", " << std::fixed << std::setprecision(10) << b << "]";

    std::ostringstream oss_config_info;
    oss_config_info << "Config: initial_step=" << config.initial_step_percent
                    << ", x_tol=" << config.x_tol << ", zero_tol=" << config.zero_tol
                    << ", min_step=" << config.min_step_percent << ", max_step=" << config.max_step_percent
                    << ", scale_limit=" << config.step_scale_limit
                    << ", plateau_thresh=" << config.plateau_mod2_threshold
                    << ", plateau_count=" << config.plateau_count_before_jump;


    INFO("Test Case: ", test_name);
    INFO(oss_interval_info.str());
    INFO(oss_config_info.str());


    finder.findRoots(a, b, found_roots);

    // --- MODIFICATION: Unconditionally print eval count to stdout ---
    std::cout << "    Test: " << test_name
              << " | Eval Count: " << finder.evalCount()
              << " | Interval: [" << std::fixed << std::setprecision(4) << a << ", " << std::setprecision(4) << b << "]"
              << " | Roots Expected: " << expected_roots.size()
              << " | Roots Found: " << found_roots.size()
              << std::endl;
    // --- END MODIFICATION ---

    // Log with Doctest's INFO for context on failure or with -s flag
    INFO("Expected roots: ", vec_to_string(expected_roots));
    INFO("Found roots   : ", vec_to_string(found_roots));
    INFO("Eval count    : ", finder.evalCount());

    CHECK(found_roots.size() == expected_roots.size());

    if (found_roots.size() == expected_roots.size()) {
        for (size_t i = 0; i < expected_roots.size(); ++i) {
            std::ostringstream oss_root_check;
            oss_root_check << "Root mismatch at index " << i << ": expected "
                           << std::fixed << std::setprecision(10) << expected_roots[i]
                           << ", found " << std::fixed << std::setprecision(10) << found_roots[i]
                           << ", difference: " << std::scientific << (found_roots[i] - expected_roots[i]);
            CHECK_MESSAGE(std::abs(found_roots[i] - expected_roots[i]) < root_comparison_tol, oss_root_check.str());
        }
    }

    if (check_eval_count) { // check_eval_count is the bool flag, max_expected_evals is the value
        CHECK(finder.evalCount() <= max_expected_evals);
    }
}



TEST_SUITE("AdaptiveBracketRootFinder Tests") {
  AdaptiveBracketConfig defaultConfig; // Default constructor values from struct definition

  // function that takes in a vector of real functions, and computes Product[1 + Exp[I fList[[i]][t]], {i, Length[fList]}]
  std::complex<double> product_exp(const std::vector<std::function<double(double)>>& fList, double t) {
    std::complex<double> result(1.0, 0.0);
    for (int i = 0; i < fList.size(); ++i) {
      result *= (1.0 + std::exp(std::complex<double>(0.0, M_PI * fList[i](t))));
    }
    return result;
  }

  AdaptiveBracketRootFinder::EvalFn f_product_exp_1 =
      [](double t) {
        std::vector<std::function<double(double)>> fList = {
          [](double t) { return t; },
          [](double t) { return t/2.0; },
      };
        return product_exp(fList, t);
  };

  TEST_CASE("Complex Product Exponential: {f(t) = t, t/2}") {
    AdaptiveBracketConfig cfg = defaultConfig;
    cfg.zero_tol = 1e-12;
    cfg.initial_step_percent = 5e-2;
    cfg.max_step_percent = 1e-1;
    check_roots("Product Exp: (1+exp(i*pi*t))", f_product_exp_1, 0.5, 3.5, {1.0, 2.0, 3.0}, cfg);
  }

  AdaptiveBracketRootFinder::EvalFn f_product_exp_2 =
      [](double t) {
        std::vector<std::function<double(double)>> fList = {
          [](double t) { return -(t-5)*(t-5)*(t-5); },
          [](double t) { return 0.01*sin(8*(t-3)); },
      };
        return product_exp(fList, t);
  };

  TEST_CASE("Complex Product Exponential w/ Plateau: {f(t) = -(t-5)^3, 0.01*sin(8*(t-3))}") {
    AdaptiveBracketConfig cfg = defaultConfig;
    cfg.zero_tol = 1e-12;
    cfg.initial_step_percent = 0.1;
    cfg.max_step_percent = 1e-1;
    cfg.min_step_percent = 1e-3;
    cfg.step_scale_limit = 5.0;
    cfg.plateau_mod2_threshold = 1e-1;
    check_roots("Product Exp w/ Platuea: (1+exp(i*pi*f(t)))",
                f_product_exp_2, 3.8, 6.2, {4.0, 6.0}, cfg);
  }

  AdaptiveBracketRootFinder::EvalFn f_product_exp_3 =
      [](double t) {
        std::vector<std::function<double(double)>> fList = {
          [](double t) { return exp(-(t-1)); },
          [](double t) { return tan(t-3); },
            [](double t) { return exp(t-5)/10; },
            [](double t) { return -exp(-(t-3)*(t-3)); }
      };
        return product_exp(fList, t);
  };

  TEST_CASE("Complex Product Exponential w/ 4 functions: {f(t) = exp(-(t-1)), tan(t-3), exp(t-5)/10, -exp(-(t-3)^2)}") {
    AdaptiveBracketConfig cfg = defaultConfig;
    cfg.zero_tol = 1e-12;
    // cfg.initial_step_percent = 5e-2;
    // cfg.max_step_percent = 1e-1;
    // there's a touching root at 3.0, but in our situation the im part crosses (at least thats what we thought, have to wait and see
    // if we run into a case where it doesn't)
    check_roots("Product Exp w/ 4 functions: (1+exp(i*pi*f(t)))",
                f_product_exp_3, 2.10, 3.90, {2.2146018366025516904, 3.7853981633974483096}, cfg);
  }
}
