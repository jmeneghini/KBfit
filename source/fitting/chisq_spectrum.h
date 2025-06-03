#ifndef CHISQ_SPECTRUM_H
#define CHISQ_SPECTRUM_H

#include "box_quant.h"
#include "chisq_base.h"
#include "task_utils.h"
#include "kbobs_handler.h"
#include "matrix.h"
#include "xml_handler.h"

class SpectrumFit : public ChiSquare {

  KBObsHandler* KBOH;
  std::vector<BoxQuantization*> BQ;
  std::vector<RVector> Ecm_over_mref;
  std::vector<RVector> m_over_mref;
  std::vector<RVector> L_mref;
  std::vector<uint> ensemble_id;
  KtildeMatrixCalculator* Kmat;
  KtildeInverseCalculator* Kinv;
  double omega_mu;

public:
  SpectrumFit(XMLHandler& xmlin, KBObsHandler* kboh,
                         XMLHandler& xmlout,
                         const std::string& outfile_stub = "");

  virtual ~SpectrumFit();

  void clear();

  void guessInitialFitParamValues(std::vector<double>& fitparams) const
  override;

  void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const override;

  void do_output(XMLHandler& xmlout) const override;

private:
  void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams)
  override;

  // Calculated once at start, then ChiSquare will use it for remainder
  // of the fit. The evalResidualsAndInvCovCholesky function just
  // updates the vars in detres; we omit that here.
  void initializeInvCov();

  friend class TaskHandler;
};

#endif //CHISQ_SPECTRUM_H
