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

  // --- Hot path data (accessed frequently in residuals calculation) ---
  // blocks - most frequently accessed
  std::vector<BoxQuantization*> BQ_blocks;//
  std::vector<uint> n_energies_per_block;//
  uint n_blocks = 0;//

  // fixed parameter flags and values (hot path)
  bool is_length_fixed = false; //
  std::vector<bool> is_mass_fixed;
  double fixed_length_value;                  // Fixed length value if is_length_fixed is true
  std::vector<double> fixed_mass_values;      // Fixed mass values indexed by decay channel*2 + particle

  // --- observations (medium frequency access) ---
  // energy data
  std::vector<RVector> Elab_samples;
  std::vector<RVector> dElab_samples;//
  std::vector<RVector> Ecm_samples;
  EnergyType residual_energy_type = EnergyType::dElab;//
  std::vector<MCObsInfo> energy_obs_infos;//

  // length and mass data
  RVector length_samples;                     // Only if not fixed (observable)
  std::vector<RVector> mass_samples;          // Only non-fixed masses (observables)

  // --- Cold path data (accessed infrequently) ---
  // prior info
  MCObsInfo length_prior;
  std::vector<MCObsInfo> mass_priors;

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
  
  std::vector<bool> are_decay_channels_identical;

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
  // Performance-critical method called for every function evaluation during fitting
  void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams)
  override;

  // Calculated once at start, then ChiSquare will use it for remainder
  // of the fit. The evalResidualsAndInvCovCholesky function just
  // updates the vars in detres; we omit that here.
  void initializeInvCovCholesky();

  friend class TaskHandler;
};

#endif //CHISQ_SPECTRUM_H
