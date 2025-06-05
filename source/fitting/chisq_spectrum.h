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
  std::vector<std::vector<RVector>> energy_samples_per_ensemble;//
  std::vector<std::vector<RVector>> mass_samples_per_ensemble;//
  std::vector<RVector> length_samples_per_ensemble;//

  AdaptiveBracketConfig root_finder_config;//
  double Elab_max, Elab_min;

  // reference length and masses for each ensemble
  // each ensemble gets a single length and multiple decay masses
  std::vector<MCObsInfo> prior_obs_infos;//
  std::vector<KBObsInfo> prior_kobs_infos;//
  std::vector<KBObsInfo> energy_kobs_infos;//

  KtildeMatrixCalculator* Kmat;//
  KtildeInverseCalculator* Kinv;//
  double omega_mu;//
  std::vector<uint> n_energies_per_block;//
  std::vector<bool> are_decay_channels_identical;//
  std::vector<uint> ensemble_id_per_block;
  std::vector<uint> n_blocks_per_ensemble;
  uint n_decay_channels;//
  uint n_kmat_params;//

public:
  SpectrumFit(XMLHandler& xmlin, KBObsHandler* kboh,
                         XMLHandler& xmlout,
                         const std::string& outfile_stub = "");

  virtual ~SpectrumFit();

  void clear();

  // order of fitinfos and fitparams are the same
  void guessInitialFitParamValues(std::vector<double>& fitparams) const
  override;

  void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const override;

  void do_output(XMLHandler& xmlout) const override;

private:
  void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams)
  override;

  void initializeFitParamsAndObservables();

  // Calculated once at start, then ChiSquare will use it for remainder
  // of the fit. The evalResidualsAndInvCovCholesky function just
  // updates the vars in detres; we omit that here.
  void initializeInvCovCholesky();

  friend class TaskHandler;
};

#endif //CHISQ_SPECTRUM_H
