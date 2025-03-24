#include "chisq_base.h"
using namespace std;

// *************************************************************

void ChiSquare::initialize_base(uint number_fit_parameters,
                                uint number_residuals,
                                uint number_resamplings) {
  if (number_fit_parameters < 1)
    throw(std::invalid_argument("Number of fit parameters > 0 required"));
  if (number_residuals <= number_fit_parameters)
    throw(std::invalid_argument(
        "Number of residuals > number fit parameters required"));
  nresiduals = number_residuals;
  nfitparams = number_fit_parameters;
  inv_cov_cholesky.resize(nresiduals);
  residuals.resize(nresiduals);
  nresamplings = number_resamplings;
}

std::string ChiSquare::output(int indent) const {
  XMLHandler xmlout;
  do_output(xmlout);
  return xmlout.output(indent);
}

std::string ChiSquare::str() const {
  XMLHandler xmlout;
  do_output(xmlout);
  return xmlout.str();
}

void ChiSquare::setResamplingIndex(uint in_resampling_index) {
  if (in_resampling_index > nresamplings)
    throw(std::invalid_argument("Invalid resampling index in ChiSquare"));
  resampling_index = in_resampling_index;
}

double ChiSquare::evalChiSquare(const vector<double>& fitparams) {
  evalResidualsAndInvCovCholesky(fitparams);
  double chisq = 0.0;
  uint nres = residuals.size();
  for (int i = nres - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    chisq += tmp * tmp;
  }
  return chisq;
}

void ChiSquare::evalDiagonalResiduals(const vector<double>& fitparams,
                                      vector<double>& diag_residuals) {
  diag_residuals.resize(nresiduals);
  evalResidualsAndInvCovCholesky(fitparams);
  uint nres = residuals.size();
  for (int i = nres - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    diag_residuals[i] = tmp;
  }
}

void ChiSquare::evalDiagonalResiduals(const vector<double>& fitparams,
                                      double* diag_residuals) {
  evalResidualsAndInvCovCholesky(fitparams);
  for (int i = nresiduals - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    diag_residuals[i] = tmp;
  }
}

// *************************************************************
