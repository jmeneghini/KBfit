/**
 * @file chisq_base.h
 * @brief Base class for correlated least-squares fitting in KBfit
 * 
 * This header defines the ChiSquare base class that provides the fundamental framework
 * for correlated χ² fitting in KBfit. All fitting methods (spectrum, determinant residual)
 * inherit from this base class.
 * 
 * Key features:
 * - Pure virtual interface for derived classes
 * - Cholesky decomposition for efficient χ² computation
 * - XML-based configuration and output
 * 
 * @author KBfit Development Team
 * @date 2024
 */

#ifndef CHISQ_BASE_H
#define CHISQ_BASE_H

#include "box_quant.h"
#include "matrix.h"
#include "mcobs_info.h"
#include "xml_handler.h"
#include "kbobs_info.h"
#include "ensemble_info.h"

/**
 * @brief Base class for correlated least-squares fitting in KBfit
 * 
 * @author KBfit Development Team
 * @date 2024
 * 
 * The ChiSquare class provides the fundamental framework for correlated χ² fitting
 * in KBfit. It cannot be constructed directly as it contains pure virtual functions
 * that must be implemented by derived classes.
 * 
 * ## Mathematical Framework
 * 
 * The fitting uses a correlated chi-square cost function:
 * 
 *     residual[j] = model[j] - obs[j]
 *     χ² = Σ_j,k residual[j] × inv_cov(j,k) × residual[k]
 * 
 * Where inv_cov is the inverse of the covariance matrix:
 * 
 *     cov(j,k) = cov(residual[j], residual[k])
 * 
 * The model can depend on the observations. A Cholesky decomposition is used:
 * 
 *     inv_cov = transpose(L) × L,   where L is lower triangular
 * 
 * ## Usage Requirements
 * 
 * Any derived class constructor must:
 * 1. Call initialize_base(number_fit_parameters, number_residuals, number_resamplings)
 *    - Sets nresiduals, nfitparams, nresamplings
 *    - Sizes residuals and inv_cov_cholesky appropriately
 * 2. Implement the required virtual functions
 * 
 * ## Required Virtual Functions
 * 
 * Derived classes must implement:
 * - do_output(XMLHandler& xmlout): Generate XML output
 * - evalResidualsAndInvCovCholesky(fitparams): Compute residuals and covariance
 * - guessInitialFitParamValues(fitparams): Provide initial parameter estimates
 * - getFitParamMCObsInfo(fitinfos): Get observable information for parameters
 * 
 * ## Performance Considerations
 * 
 * The evalResidualsAndInvCovCholesky method is called repeatedly during fitting
 * and should be optimized for performance. The base class provides efficient
 * χ² computation using the pre-computed Cholesky decomposition.
 */

class ChiSquare {

private:
#ifndef NO_CXX11
  ChiSquare(const ChiSquare&) = delete;              ///< Copy constructor disabled
  ChiSquare& operator=(const ChiSquare&) = delete;   ///< Assignment operator disabled
#else
  ChiSquare(const ChiSquare&);                       ///< Copy constructor disabled (C++03)
  ChiSquare& operator=(const ChiSquare&);            ///< Assignment operator disabled (C++03)
#endif

protected: 
  /// @name Core Fitting Data
  /// @{
  uint nresiduals;                                   ///< Number of residuals in the fit
  uint nfitparams;                                   ///< Number of fit parameters
  LowerTriangularMatrix<double> inv_cov_cholesky;    ///< Cholesky decomposition of inverse covariance
  std::vector<double> residuals;                     ///< Current residual values
  uint nresamplings;                                 ///< Number of bootstrap/jackknife resamplings
  uint resampling_index;                             ///< Current resampling index (0 = full sample)
  BoxQuantization::QuantCondType qctype_enum;        ///< Type of quantization condition
  /// @}

public:
  /// @name Construction and Destruction
  /// @{
  ChiSquare() {}                                     ///< Default constructor
  virtual ~ChiSquare() {}                            ///< Virtual destructor
  /// @}

  /// @name Accessors
  /// @{
  uint getNumberOfFitParameters() const { return nfitparams; }   ///< Get number of fit parameters
  uint getNumberOfResiduals() const { return nresiduals; }       ///< Get number of residuals
  uint getNumberOfResamplings() const { return nresamplings; }   ///< Get number of resamplings
  /// @}

  /// @name Resampling Control
  /// @{
  /**
   * @brief Set the current resampling index
   * @param in_resampling_index Index for bootstrap/jackknife (0 = full sample)
   */
  void setResamplingIndex(uint in_resampling_index);
  /// @}

  /// @name Output Functions
  /// @{
  /**
   * @brief Generate XML output (pure virtual - must be implemented by derived classes)
   * @param xmlout XML handler for output
   */
  virtual void do_output(XMLHandler& xmlout) const = 0;

  /**
   * @brief Generate formatted XML output string
   * @param indent Indentation level
   * @return Formatted XML string
   */
  std::string output(int indent = 0) const;

  /**
   * @brief Generate compact XML output string
   * @return Compact XML string
   */
  std::string str() const;
  /// @}

  /// @name Pure Virtual Interface
  /// @{
  /**
   * @brief Generate initial parameter guesses (pure virtual)
   * @param fitparams Vector to store initial parameter values
   * @param only_update_priors If true, only update prior-related parameters
   */
  virtual void guessInitialFitParamValues(std::vector<double>& fitparams, bool only_update_priors) const = 0;

  /**
   * @brief Get observable information for fit parameters (pure virtual)
   * @param fitinfos Vector to store MCObsInfo for each parameter
   */
  virtual void getFitParamMCObsInfo(std::vector<MCObsInfo>& fitinfos) const = 0;
  /// @}

  /// @name Chi-Square Evaluation
  /// @{
  /**
   * @brief Evaluate chi-square value for given parameters
   * @param fitparams Current fit parameter values
   * @return Chi-square value
   */
  double evalChiSquare(const std::vector<double>& fitparams);

  /**
   * @brief Evaluate diagonal residuals (whitened residuals)
   * @param fitparams Current fit parameter values
   * @param diag_residuals Vector to store diagonal residuals
   */
  void evalDiagonalResiduals(const std::vector<double>& fitparams,
                             std::vector<double>& diag_residuals);

  /**
   * @brief Evaluate diagonal residuals (whitened residuals)
   * @param fitparams Current fit parameter values  
   * @param diag_residuals Array to store diagonal residuals
   */
  void evalDiagonalResiduals(const std::vector<double>& fitparams,
                             double* diag_residuals);
  /// @}

  /// @name Static Utility Functions
  /// @{
  /**
   * @brief Read observable information from XML with optional fixed values
   * @param xmlin XML handler for input
   * @param tag XML tag to search within
   * @param get_name Whether to read name attribute
   * @param obskey Output observable key
   * @param kset Set of observable keys (for duplicate checking)
   * @param name Output name string
   * @param mcens Monte Carlo ensemble information
   * @param fixed_values Map of fixed parameter values
   */
  static void read_obs(XMLHandler& xmlin, const std::string& tag, bool get_name,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset,
                       std::string& name, const MCEnsembleInfo& mcens,
                       std::map<KBObsInfo, double>& fixed_values);

  /**
   * @brief Read observable information from XML (simple version)
   * @param xmlin XML handler for input
   * @param tag XML tag to search within
   * @param obskey Output observable key
   * @param kset Set of observable keys (for duplicate checking)
   */
  static void read_obs(XMLHandler& xmlin, const std::string& tag,
                       MCObsInfo& obskey, std::set<MCObsInfo>& kset);
  /// @}

protected:
  /// @name Protected Interface for Derived Classes
  /// @{
  /**
   * @brief Initialize base class data structures
   * @param number_fit_parameters Number of fit parameters
   * @param number_residuals Number of residuals
   * @param number_resamplings Number of bootstrap/jackknife resamplings
   * 
   * This method must be called by derived class constructors to set up
   * the base class data structures properly.
   */
  void initialize_base(uint number_fit_parameters, uint number_residuals,
                       uint number_resamplings);

  /**
   * @brief Evaluate residuals and inverse covariance Cholesky decomposition (pure virtual)
   * @param fitparams Current fit parameter values
   * 
   * This is the core method that derived classes must implement. It should:
   * 1. Compute model predictions for the given parameters
   * 2. Calculate residuals = predictions - observations
   * 3. Store residuals in the residuals member variable
   * 
   * The covariance matrix and its Cholesky decomposition are typically
   * computed once during initialization and stored in inv_cov_cholesky.
   * 
   * @note This method is called repeatedly during fitting and should be optimized.
   */
  virtual void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) = 0;
  /// @}


};

// *****************************************************************
#endif
