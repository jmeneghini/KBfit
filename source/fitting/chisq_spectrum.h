#ifndef CHISQ_SPECTRUM_H
#define CHISQ_SPECTRUM_H

#include "box_quant.h"
#include "chisq_base.h"
#include "kbobs_handler.h"
#include "matrix.h"
#include "xml_handler.h"

enum class EnergyType {
  Elab,
  Ecm,
  dElab
};

struct EnsembleFitData {
  // this ensemble_info and id - default to independent
  MCEnsembleInfo ensemble_info = MCEnsembleInfo(0);

  // blocks
  std::vector<BoxQuantization*> BQ_blocks = {}; // Box quantization blocks
  std::vector<uint> n_energies_per_block = {}; // Number of energies per block
  std::vector<std::pair<double, double>> Ecm_bounds_per_block = {};
  uint n_blocks = 0;//

  // lattice parameters
  uint Llat = 0; // Lattice spatial extent

  // fixed parameter flags and values
  bool is_anisotropy_fixed = true; // Default to isotropic (fixed anisotropy = 1.0)
  std::vector<bool> is_mass_fixed = {}; // Indexed by decay channel*2 + particle
  double fixed_anisotropy_value = 1.0;         // Fixed anisotropy value (default 1.0 for isotropic)
  std::vector<double> fixed_mass_values = {};  // Fixed mass values indexed by decay channel*2 + particle

  // energy data
  std::vector<RVector> Elab_samples = {}; // Elab samples for each ensemble
  std::vector<RVector> dElab_samples = {};
  std::vector<RVector> Ecm_samples = {};
  EnergyType residual_energy_type = EnergyType::dElab;//
  std::vector<MCObsInfo> energy_obs_infos = {};
  std::vector<NonInteractingPair> non_interacting_pairs = {}; // Non-interacting pairs for each energy

  // reference mass and mass data (mref is always a parameter)
  RVector mref_samples;                       // Always present (KBScale observable)
  RVector anisotropy_samples;                 // Only if not fixed (anisotropic case)
  std::vector<RVector> mass_samples = {};     // Only non-fixed masses (observables)

  // prior info
  MCObsInfo mref_prior = MCObsInfo("KBScale", 0);
  MCObsInfo anisotropy_prior = MCObsInfo("default", 0);
  std::vector<MCObsInfo> mass_priors = {}; // Indexed by decay channel*2 + particle
};


class SpectrumFit : public ChiSquare {

  KBObsHandler* KBOH;

  AdaptiveBracketConfig root_finder_config;//

  KtildeMatrixCalculator* Kmat;//
  KtildeInverseCalculator* Kinv;//
  uint n_kmat_params;//
  uint n_decay_channels;

  double omega_mu;//

  std::vector<EnsembleFitData> ensemble_fit_data;
  
  std::vector<bool> are_decay_channels_identical;

  // Reusable temporary vectors for performance (pre-allocated)
  mutable std::vector<double> energy_shift_predictions;
  mutable std::vector<std::pair<double, NonInteractingPair>> shift_obs_w_NIs;
  mutable std::vector<uint> fn_calls;
  mutable std::vector<std::pair<double, double>> decay_channel_masses;

public:
  SpectrumFit(XMLHandler& xmlin, KBObsHandler* kboh,
                         XMLHandler& xmlout,
                         const std::string& outfile_stub = "");

  virtual ~SpectrumFit();

  void clear();

  // order of fitinfos and fitparams are the same
  void guessInitialFitParamValues(std::vector<double>& fitparams, bool only_update_priors) const
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

  // Helper function to parse non-interacting pair strings
  NonInteractingPair parseNonInteractingPair(const std::string& pair_str, 
                                             const std::vector<DecayChannelInfo>& decay_channels) const;

  friend class TaskHandler;
};

#endif //CHISQ_SPECTRUM_H
