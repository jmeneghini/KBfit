#ifndef CHISQ_SPECTRUM_H
#define CHISQ_SPECTRUM_H

/**
 * @file chisq_spectrum.h
 * @brief Spectrum fitting implementation for K-matrix parameter estimation in lattice QCD
 * 
 * This file contains the SpectrumFit class which implements spectrum fitting methodology
 * for two-hadron systems in finite volume using the Lüscher method. The spectrum fitting
 * approach solves the quantization condition QC(Ecm) = 0 to find energy eigenvalues,
 * then minimizes χ² residuals between model predictions and observed energies.
 * 
 * Key features:
 * - Root finding for quantization conditions (always uses StildeCB)
 * - χ² minimization with correlated error analysis
 * - Support for multiple ensembles and momentum blocks
 * - MPI parallelization for resampling fits
 * - Performance optimizations for memory access and computational efficiency
 * 
 * @author KBfit Development Team
 * @date 2024
 */

#include "box_quant.h"
#include "chisq_base.h"
#include "kbobs_handler.h"
#include "matrix.h"
#include "xml_handler.h"
#include <memory>

/**
 * @enum EnergyType
 * @brief Enumeration for different energy representations in spectrum fitting
 * 
 * Defines the various energy types used in two-hadron scattering analysis:
 * - Elab: Laboratory frame energy
 * - Ecm: Center-of-mass frame energy  
 * - dElab: Energy shift relative to non-interacting reference
 */
enum class EnergyType {
  Elab,   ///< Laboratory frame energy
  Ecm,    ///< Center-of-mass frame energy
  dElab   ///< Energy shift (difference from non-interacting)
};

/**
 * @struct EnsembleFitData
 * @brief Container for ensemble-specific data required for spectrum fitting
 * 
 * This structure holds all the necessary data for performing spectrum fits on a single
 * ensemble, including energy levels, lattice parameters, mass information, and prior
 * constraints. The data is organized to support efficient access during the intensive
 * fitting process.
 */
struct EnsembleFitData {
  /// Ensemble information and metadata (defaults to independent ensemble)
  MCEnsembleInfo ensemble_info = MCEnsembleInfo(0);

  /// @name Box Quantization Data
  /// @{
  std::vector<BoxQuantization*> BQ_blocks = {};              ///< Box quantization blocks for different momentum classes
  std::vector<uint> n_energies_per_block = {};              ///< Number of energy levels per momentum block
  std::vector<std::pair<double, double>> Ecm_bounds_per_block = {}; ///< Energy bounds for each block
  uint n_blocks = 0;                                         ///< Total number of momentum blocks
  /// @}

  /// @name Lattice Parameters
  /// @{
  uint Llat = 0;                                            ///< Lattice spatial extent
  /// @}

  /// @name Fixed Parameters Configuration
  /// @{
  bool is_anisotropy_fixed = true;                          ///< Whether anisotropy is fixed (default: isotropic)
  std::vector<bool> is_mass_fixed = {};                     ///< Fixed mass flags (indexed by decay channel*2 + particle)
  double fixed_anisotropy_value = 1.0;                      ///< Fixed anisotropy value (default: 1.0 for isotropic)
  std::vector<double> fixed_mass_values = {};               ///< Fixed mass values (indexed by decay channel*2 + particle)
  /// @}

  /// @name Energy Data
  /// @{
  std::vector<RVector> Elab_samples = {};                   ///< Laboratory frame energy samples
  std::vector<RVector> dElab_samples = {};                  ///< Energy shift samples
  std::vector<RVector> Ecm_samples = {};                    ///< Center-of-mass frame energy samples
  EnergyType residual_energy_type = EnergyType::dElab;      ///< Type of energy used for residual calculation
  std::vector<MCObsInfo> energy_obs_infos = {};             ///< Observable information for each energy level
  std::vector<NonInteractingPair> non_interacting_pairs = {}; ///< Non-interacting pairs for each energy
  /// @}

  /// @name Mass and Reference Data
  /// @{
  RVector mref_samples;                                     ///< Reference mass samples (always present as KBScale)
  RVector anisotropy_samples;                               ///< Anisotropy samples (only if not fixed)
  std::vector<RVector> mass_samples = {};                   ///< Mass samples for non-fixed masses
  /// @}

  /// @name Prior Information
  /// @{
  MCObsInfo mref_prior = MCObsInfo("KBScale", 0);           ///< Reference mass prior
  MCObsInfo anisotropy_prior = MCObsInfo("default", 0);     ///< Anisotropy prior
  std::vector<MCObsInfo> mass_priors = {};                  ///< Mass priors (indexed by decay channel*2 + particle)
  /// @}
};


/**
 * @class SpectrumFit
 * @brief Spectrum fitting class for K-matrix parameter estimation using the Lüscher method
 * 
 * @author KBfit Development Team
 * @date 2024
 * 
 * This class implements the spectrum fitting methodology for two-hadron systems in finite
 * volume. It differs from determinant residual (detres) fitting by solving the Omega
 * function Ω(Ecm) = 0 to predict energy shifts from non-interacting (NI) levels, then
 * performing χ² minimization on the residuals between predicted and observed energy shifts.
 * 
 * Key features:
 * - Uses StildeCB quantization condition exclusively
 * - Root finding to solve Omega function for energy shift predictions
 * - Supports multiple ensembles and momentum blocks (AR, OA, PD, CD)
 * - Incorporates priors for non-dElab observables (mass, anisotropy)
 * - Thread-safe design for MPI parallelization
 * 
 * @see ChiSquare Base class for correlated χ² fitting
 * @see BoxQuantization For quantization condition calculations
 * @see KtildeMatrixCalculator For K-matrix computations
 */
class SpectrumFit : public ChiSquare {
private:
  /// @name Core Components
  /// @{
  KBObsHandler* KBOH;                                       ///< Observable handler for data management
  AdaptiveBracketConfig root_finder_config;                ///< Configuration for root finding algorithms
  /// @}

  /// @name K-matrix Calculation
  /// @{
  KtildeMatrixCalculator* Kmat;                            ///< K-matrix calculator
  KtildeInverseCalculator* Kinv;                           ///< K-matrix inverse calculator
  uint n_kmat_params;                                      ///< Number of K-matrix parameters
  uint n_decay_channels;                                   ///< Number of decay channels
  /// @}

  /// @name Physics Parameters
  /// @{
  double omega_mu;                                         ///< Omega parameter for scattering calculations
  /// @}

  /// @name Ensemble Data
  /// @{
  std::vector<EnsembleFitData> ensemble_fit_data;          ///< Data for all ensembles in the fit
  std::vector<bool> are_decay_channels_identical;          ///< Flag for identical decay channels
  /// @}

  /// @name Performance Optimization
  /// @{
  /// Pre-allocated temporary vectors to avoid repeated memory allocations during fitting.
  /// These are marked mutable to allow modification in const methods during fitting.
  mutable std::vector<double> energy_shift_predictions;                    ///< Cached energy shift predictions
  mutable std::vector<std::pair<double, NonInteractingPair>> shift_obs_w_NIs; ///< Cached observable/NI pairs
  mutable std::vector<uint> fn_calls;                                      ///< Cached function call counts
  mutable std::vector<std::pair<double, double>> decay_channel_masses;     ///< Cached decay channel masses
  /// @}

public:
  /// @name Construction and Destruction
  /// @{
  
  /**
   * @brief Constructor for spectrum fitting from XML configuration
   * @param xmlin XML handler for input configuration
   * @param kboh Observable handler for data management
   * @param xmlout XML handler for output
   * @param outfile_stub Output file stub for result files
   */
  SpectrumFit(XMLHandler& xmlin, KBObsHandler* kboh,
              XMLHandler& xmlout,
              const std::string& outfile_stub = "");

  /**
   * @brief Virtual destructor
   */
  virtual ~SpectrumFit();

  /**
   * @brief Clear all internal data structures
   */
  void clear();

  /**
   * @brief Deep copy/clone method to create an identical object with new pointers
   * @param new_kboh New observable handler (optional, uses existing if nullptr)
   * @return Unique pointer to cloned SpectrumFit object
   */
  std::unique_ptr<SpectrumFit> clone(KBObsHandler* new_kboh = nullptr) const;
  /// @}

  /// @name ChiSquare Interface Implementation
  /// @{
  
  /**
   * @brief Generate initial guess for fit parameters
   * @param fitparams Vector to store initial parameter values
   * @param only_update_priors If true, only update prior-related parameters
   * 
   * Order of fitinfos and fitparams are consistent between methods.
   */
  void guessInitialFitParamValues(std::vector<double>& fitparams, bool only_update_priors) const override;

  /**
   * @brief Get observable information for all fit parameters
   * @param fitinfos Vector to store MCObsInfo for each parameter
   */
  void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const override;

  /**
   * @brief Generate XML output for fit results
   * @param xmlout XML handler for output
   */
  void do_output(XMLHandler& xmlout) const override;
  /// @}

private:
  /// @name Core Fitting Methods
  /// @{
  
  /**
   * @brief Performance-critical method called for every function evaluation during fitting
   * @param fitparams Current fit parameter values
   * 
   * This method computes residuals between predicted and observed energy shifts from
   * non-interacting levels. It solves the Omega function Ω(Ecm) = 0 via root
   * finding to predict energy shifts, then calculates χ² residuals.
   */
  void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) override;

  /**
   * @brief Private default constructor for clone method
   */
  SpectrumFit();

  /**
   * @brief Initialize inverse covariance Cholesky decomposition
   * 
   * Calculated once at start, then used by ChiSquare base class for χ² computation.
   */
  void initializeInvCovCholesky();

  /**
   * @brief Helper function to parse non-interacting pair strings
   * @param pair_str String representation of non-interacting pair
   * @param decay_channels Vector of decay channel information
   * @return Parsed NonInteractingPair object
   */
  NonInteractingPair parseNonInteractingPair(const std::string& pair_str, 
                                             const std::vector<DecayChannelInfo>& decay_channels) const;
  /// @}

  friend class TaskHandler;
};

#endif //CHISQ_SPECTRUM_H
