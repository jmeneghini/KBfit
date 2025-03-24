#include "minimizer.h"
#include "xml_handler.h"
using namespace std;

class DummyChiSquare : public ChiSquare {
  vector<double> exact_answer;
  vector<double> initial_guess;

public:
  DummyChiSquare(vector<double>& errors) : exact_answer(4), initial_guess(4) {
    uint numfitparams = 4;
    uint numresiduals = 12;
    uint numresamplings = 100;
    if (errors.size() != 12)
      throw(std::invalid_argument("Should give 12 errors"));
    exact_answer[0] = 1.0;
    initial_guess[0] = 0.5;
    exact_answer[1] = 2.1;
    initial_guess[1] = 3.2;
    exact_answer[2] = 3.0;
    initial_guess[2] = 1.4;
    exact_answer[3] = 0.5;
    initial_guess[3] = 0.4;
    initialize_base(numfitparams, numresiduals, numresamplings);
    for (uint i = 0; i < nresiduals; i++) {
      inv_cov_cholesky(i, i) = 1.0 / (errors.at(i) * errors.at(i));
      for (uint j = 0; j < i; ++j)
        inv_cov_cholesky(i, j) = 0.0;
    }
  }

  virtual ~DummyChiSquare() {}

  /*  double model(const vector<double>& params, double t) const
     { return params[0]*exp(-params[1]*double(t))*(1.0
             +params[2]*exp(-params[3]*params[3]*double(t)));} */

  double model(const vector<double>& params, double t) const {
    return (t - params[0]) * (t - params[1]) * (t - params[2]) * params[3];
  }

  virtual void evalResidualsAndInvCovCholesky(const vector<double>& fitparams) {
    double chisq = 0.0;
    for (int t = 0; t < int(nresiduals); t++) {
      residuals[t] = model(fitparams, 0.5 * double(t)) -
                     model(exact_answer, 0.5 * double(t));
      chisq += residuals[t] * residuals[t];
    }
    cout << "chisq=" << chisq << endl;
  }

  virtual void guessInitialFitParamValues(vector<double>& fitparams) const {
    fitparams = initial_guess;
  }

  virtual void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfo) const {}

  virtual void do_output(XMLHandler& xmlout) const {}
};

void testMinimizer(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMinimizer") == 0)
    return;

  cout << endl << "Starting TestMinimizer" << endl;

  XMLHandler xmlr(xml_in, "TestMinimizer");

  ChiSquareMinimizerInfo csqinfo(xmlr);
  vector<double> errors;
  xmlread(xmlr, "Errors", errors, "TestMinimizer");
  DummyChiSquare dummy(errors);

  ChiSquareMinimizer CSM(dummy, csqinfo);

  vector<double> init_params;
  dummy.guessInitialFitParamValues(init_params);
  double chisq;
  vector<double> params_at_min;
  RealSymmetricMatrix params_covariance;
  XMLHandler xmlout;

  bool flag = CSM.findMinimum(init_params, chisq, params_at_min,
                              params_covariance, xmlout);
  cout << "flag = " << flag << endl;
  cout << xmlout.output() << endl;
  for (uint k = 0; k < params_at_min.size(); ++k)
    cout << "params_at_min[" << k << "] = " << params_at_min[k] << endl;
  for (uint k = 0; k < params_covariance.size(); ++k)
    for (uint j = k; j < params_covariance.size(); ++j)
      cout << "params_cov[" << k << "," << j
           << "] = " << params_covariance(k, j) << endl;
  cout << "chisq = " << chisq << endl;
}
