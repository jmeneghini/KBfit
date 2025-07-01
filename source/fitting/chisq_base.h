#ifndef CHISQ_BASE_H
#define CHISQ_BASE_H

#include "box_quant.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "xml_handler.h"
#include "kbobs_info.h"
#include "ensemble_info.h"

// ********************************************************************************
// * *
// *   "ChiSquare" is the all-important base class for correlated least-squares
// *
// *   fitting in KBfit.  The base class cannot be constructed on its own since
// *
// *   it contains several purely virtual member functions. *
// * *
// *   For the fitting, the "cost" function is a correlated chi-square: *
// * *
// *      residual[j] = model[j]-obs[j] *
// * *
// *      chi-square = sum_j residual[j] * inv_cov(j,k) * residual[k] *
// * *
// *   where inv_cov = the inverse of the covariance matrix given by *
// * *
// *      cov(j,k) = cov( residual[j], residual[k] ) *
// * *
// *   The model could depend on the observations. *
// * *
// *   A Cholesky decomposition of inv_cov can be done: *
// *              inv_cov = transpose(L) * L,   L = lower triangular *
// *   The fit parameters determine the model[j] values, and the fit parameters
// *
// *   are adjusted until the chi-square is a minimum. *
// * *
// * *
// *   The constructor of any class "DerivedFit" derived from "ChiSquare" must *
// * *
// *     --  call "initialize_base" *
// *           initialize_base(number_fit_parameters, number_residuals, *
// *                           number_resamplings); *
// *       which -- sets "nresiduals" the number of residuals *
// *             -- sets "nfitparams" the number of fit parameters *
// *             -- sets "nresamplings" the number of resamplings *
// *             -- appropriately sizes "residuals" and "inv_cov_cholesky" *
// * *
// *   The derived class must define the members below. *
// * *
// *        virtual void output(XMLHandler& xmlout) const; *
// * *
// *        virtual void evalResidualsAndInvCovCholesky( *
// *                           const std::vector<double>& fitparams); *
// * *
// *        virtual void guessInitialFitParamValues( *
// *                           std::vector<double>& fitparams) const; *
// * *
// * *
// ********************************************************************************

class ChiSquare {

private:
#ifndef NO_CXX11
  ChiSquare(const ChiSquare&) = delete;
  ChiSquare& operator=(const ChiSquare&) = delete;
#else
  ChiSquare(const ChiSquare&);
  ChiSquare& operator=(const ChiSquare&);
#endif

protected: // derived classes have access to the protected members
  uint nresiduals;
  uint nfitparams;
  LowerTriangularMatrix<double> inv_cov_cholesky;
  std::vector<double> residuals;
  uint nresamplings;
  uint resampling_index;
  BoxQuantization::QuantCondType qctype_enum;

public:
  ChiSquare() {}

  virtual ~ChiSquare() {}

  uint getNumberOfFitParameters() const { return nfitparams; }

  uint getNumberOfResiduals() const { return nresiduals; }

  uint getNumberOfResamplings() const { return nresamplings; }

  void setResamplingIndex(uint in_resampling_index);

  virtual void do_output(XMLHandler& xmlout) const = 0;

  std::string output(int indent = 0) const;

  std::string str() const;

  virtual void
  guessInitialFitParamValues(std::vector<double>& fitparams, bool only_update_priors) const = 0;

  virtual void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const = 0;

  double evalChiSquare(const std::vector<double>& fitparams);

  void evalDiagonalResiduals(const std::vector<double>& fitparams,
                             std::vector<double>& diag_residuals);

  void evalDiagonalResiduals(const std::vector<double>& fitparams,
                             double* diag_residuals);

  static void read_obs(XMLHandler& xmlin, const std::string& tag, bool get_name,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset,
                       std::string& name, const MCEnsembleInfo& mcens,
                       std::map<KBObsInfo, double>& fixed_values);

  static void read_obs(XMLHandler& xmlin, const std::string& tag,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset);


protected:
  void initialize_base(uint number_fit_parameters, uint number_residuals,
                       uint number_resamplings);

  virtual void
  evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) = 0;


};

// *****************************************************************
#endif
