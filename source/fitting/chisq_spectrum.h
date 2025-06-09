#ifndef CHISQ_SPECTRUM_H
#define CHISQ_SPECTRUM_H

#include "box_quant.h"
#include "chisq_base.h"
#include "task_utils.h"
#include "kbobs_handler.h"
#include "matrix.h"
#include "xml_handler.h"

enum class EnergyType {
  Elab,
  Ecm,
  dElab
};

struct EnsembleFitData {
  // this ensemble_info and id
  MCEnsembleInfo ensemble_info;

  // --- observations ---
  // energy data
  std::vector<RVector> Elab_samples;
  std::vector<RVector> dElab_samples;//
  std::vector<RVector> Ecm_samples;
  EnergyType residual_energy_type = EnergyType::dElab;//
  std::vector<MCObsInfo> energy_obs_infos;//

  // length data
  RVector length_samples;

  // mass data
  std::vector<RVector> mass_samples;

  // blocks
  std::vector<BoxQuantization*> BQ_blocks;//

  // prior info
  // each ensemble gets a single length and multiple decay masses
  MCObsInfo length_prior;
  bool is_length_fixed; //

  std::vector<MCObsInfo> mass_priors;
  std::vector<bool> is_mass_fixed;
  std::vector<bool> are_decay_channels_identical;

  // useful counters
  uint n_blocks;//
  std::vector<uint> n_energies_per_block;//

  // Constructor
  explicit EnsembleFitData(const MCEnsembleInfo& info)
    : ensemble_info(info) {}
};


class SpectrumFit : public ChiSquare {

  KBObsHandler* KBOH;

  AdaptiveBracketConfig root_finder_config;//
  double Elab_max, Elab_min;

  KtildeMatrixCalculator* Kmat;//
  KtildeInverseCalculator* Kinv;//
  uint n_kmat_params;//
  uint n_decay_channels;

  double omega_mu;//

  std::vector<EnsembleFitData> ensemble_fit_data;

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
