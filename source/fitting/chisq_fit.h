/**
 * @file chisq_fit.h
 * @brief Chi-square fitting workflow interface with MPI support
 *
 * This header provides the main interface for chi-square fitting in KBfit,
 * supporting both serial and MPI parallel execution modes. The fitting process
 * uses bootstrap/jackknife resampling to estimate parameter uncertainties.
 *
 * Key features:
 * - Automatic MPI detection and load balancing
 * - Progress tracking with visual progress bars
 * - Comprehensive error handling and logging
 * - Chi-square quality assessment utilities
 *
 * @author KBfit Development Team
 * @date 2024
 */

#ifndef CHISQ_FIT_H
#define CHISQ_FIT_H

#include "chisq_base.h"
#include "indicators.h" // single-header Indicators library
#include "kbobs_handler.h"
#include "mc_estimate.h"
#include "minimizer.h"
#include "param_registry.h"
#include "xml_handler.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <string>

/**
 * @brief Display progress bar for fitting operations
 * @param i Current iteration number
 * @param total Total number of iterations
 * @param is_slurm True if running in SLURM environment, false otherwise
 *
 * Shows a visual progress bar with elapsed time and estimated time remaining.
 * Used during bootstrap/jackknife resampling to provide user feedback.
 */
void show_progress(std::size_t i, std::size_t total, bool is_slurm = false);

/// @name Main Fitting Functions
/// @{

/**
 * @brief Main entry point for chi-square fitting with automatic MPI detection
 * @param chisq_ref Reference to ChiSquare object (SpectrumFit or
 * DeterminantResidualFit)
 * @param csm_info Minimizer configuration information
 * @param chisq_dof Output chi-square per degree of freedom
 * @param fitqual Output fit quality (p-value)
 * @param bestfit_params Output best-fit parameter estimates with uncertainties
 * @param param_covariance Output parameter covariance matrix
 * @param out_sampling_file Output file path for sampling results
 * @param xmlout XML handler for logging and output
 * @param kobs Observable handler for data management
 *
 * Automatically detects MPI environment and chooses appropriate execution mode:
 * - Single process: Uses serial implementation
 * - Multiple processes: Uses MPI parallel implementation with load balancing
 *
 * The fitting process:
 * 1. Performs minimization on full sample to get best-fit parameters
 * 2. Loops over bootstrap/jackknife resamplings to compute uncertainties
 * 3. Generates comprehensive output including fit quality assessment
 *
 * @throws std::invalid_argument if fitting fails or insufficient data
 */
void doChiSquareFitting(ChiSquare& chisq_ref,
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual,
                        std::vector<MCEstimate>& bestfit_params,
                        RealSymmetricMatrix& param_covariance,
                        const std::string& out_sampling_file,
                        XMLHandler& xmlout, KBObsHandler* kobs);

/**
 * @brief Serial implementation of chi-square fitting
 * @param chisq_ref Reference to ChiSquare object
 * @param csm_info Minimizer configuration
 * @param chisq_dof Output chi-square per degree of freedom
 * @param fitqual Output fit quality
 * @param bestfit_params Output best-fit parameters
 * @param param_covariance Output parameter covariance matrix
 * @param out_sampling_file Output file for sampling results
 * @param xmlout XML handler for output
 * @param kobs Observable handler
 *
 * Performs fitting in serial mode with progress tracking.
 */
void doChiSquareFittingSerial(ChiSquare& chisq_ref,
                              const ChiSquareMinimizerInfo& csm_info,
                              double& chisq_dof, double& fitqual,
                              std::vector<MCEstimate>& bestfit_params,
                              RealSymmetricMatrix& param_covariance,
                              const std::string& out_sampling_file,
                              XMLHandler& xmlout, KBObsHandler* kobs);

/**
 * @brief MPI parallel implementation of chi-square fitting
 * @param chisq_ref Reference to ChiSquare object
 * @param csm_info Minimizer configuration
 * @param chisq_dof Output chi-square per degree of freedom
 * @param fitqual Output fit quality
 * @param bestfit_params Output best-fit parameters
 * @param param_covariance Output parameter covariance matrix
 * @param out_sampling_file Output file for sampling results
 * @param xmlout XML handler for output
 * @param kobs Observable handler
 * @param comm MPI communicator
 *
 * Performs fitting with MPI parallelization:
 * - Rank 0 handles full sample fit and result aggregation
 * - All ranks participate in resampling fits with load balancing
 * - Efficient result gathering and output generation
 */
void doChiSquareFittingMPI(
    ChiSquare& chisq_ref, const ChiSquareMinimizerInfo& csm_info,
    double& chisq_dof, double& fitqual, std::vector<MCEstimate>& bestfit_params,
    RealSymmetricMatrix& param_covariance, const std::string& out_sampling_file,
    XMLHandler& xmlout, KBObsHandler* kobs, MPI_Comm comm);

/// @}

/// @name Statistical Functions
/// @{

/**
 * @brief Upper incomplete gamma function ratio
 * @param s Shape parameter (must be > 0)
 * @param x Value parameter (must be >= 0)
 * @return Q(s,x) = Γ(s,x)/Γ(s) = 1 - P(s,x)
 *
 * Q(s,0) = 1, Q(s,∞) = 0, and decreases monotonically.
 * Used for chi-square p-value calculations.
 *
 * @throws std::invalid_argument for invalid parameters
 */
double Qgamma(double s, double x);

/**
 * @brief Lower incomplete gamma function ratio
 * @param s Shape parameter (must be > 0)
 * @param x Value parameter (must be >= 0)
 * @return P(s,x) = γ(s,x)/Γ(s) = 1 - Q(s,x)
 *
 * P(s,0) = 0, P(s,∞) = 1, and increases monotonically.
 *
 * @throws std::invalid_argument for invalid parameters
 */
double Pgamma(double s, double x);

/**
 * @brief Calculate chi-square fit quality (p-value)
 * @param dof Degrees of freedom
 * @param chisquare Chi-square value
 * @return Fit quality = Q(dof/2, chisquare/2)
 *
 * Returns the probability that a chi-square random variable with
 * dof degrees of freedom would exceed the observed value.
 *
 * Interpretation:
 * - Values near 1: Very good fit
 * - Values near 0: Very poor fit
 * - Values < 0.05: Typically considered poor fit
 */
double getChiSquareFitQuality(unsigned int dof, double chisquare);

/// @}

// ****************************************************************
#endif
