/**
 * @file chisq_fit.cc
 * @brief Implementation of chi-square fitting workflow with MPI support
 *
 * This file implements the complete chi-square fitting workflow for KBfit,
 * including both serial and MPI parallel execution modes. The fitting process
 * uses bootstrap/jackknife resampling to estimate parameter uncertainties.
 *
 * Key features:
 * - Automatic detection of MPI environment
 * - Load-balanced parallel resampling fitting across MPI ranks
 * - Progress tracking with visual progress bars
 * - Comprehensive error handling and logging
 * - Efficient result aggregation and output generation
 *
 * @author KBfit Development Team
 * @date 2024
 */

#include "chisq_fit.h"
#include <algorithm>
#include <cstdio>
#include <mpi.h>
#include <sstream>
#include <thread>
#include <iostream>
#include <iomanip>
using namespace std;

/**
 * @brief Display progress bar for fitting operations
 * @param i Current iteration number
 * @param total Total number of iterations
 * @param is_slurm True if running in SLURM environment, false otherwise
 *
 * Implementation of visual progress bar with elapsed time and estimated time
 * remaining. Uses static state to persist across calls for the same fitting
 * session.
 */
void show_progress(std::size_t i, std::size_t total, bool is_slurm) {
  // Use a static pointer to persist the bar across calls.
  static std::unique_ptr<indicators::ProgressBar> bar_ptr;
  // Use a static total to detect when a new fit with a different number of
  // samples starts.
  static std::size_t bar_total = 0;
  // Track last update for SLURM buffering
  static std::size_t last_update = 0;

  // SLURM mode: Simple text-based progress updates
  if (is_slurm) {
    // Initialize on first call or total change
    if (bar_total != total) {
      bar_total = total;
      last_update = 0;
    }
    
    // Print progress every 5 percent, at start, and when complete.
    bool should_print = (i == 0) || (i >= total) ||
                            ((i - last_update) >= (total / 20)) ||
                            ((i % (total / 20)) == 0);
    
    if (should_print) {
      if (i >= total) {
        std::cerr << "Sample " << total << "/" << total << " - Fit completed successfully!" << std::endl;
      } else {
        double percent = (100.0 * i) / total;
        std::cerr << "Sample " << i << "/" << total << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::endl;
      }
      std::cerr.flush();
      last_update = i;
    }
    return;
  }
  
  // Non-SLURM mode: Use visual progress bar
  // (Re)initialize if this is the first call, or if the total number of samples has changed.
  if (!bar_ptr || bar_total != total) {
    bar_total = total;
    last_update = 0;
    
    std::cerr << "Starting fit with " << total << " samples..." << std::endl;
    std::cerr.flush();
    
    bar_ptr = std::make_unique<indicators::ProgressBar>(
        indicators::option::BarWidth{50}, indicators::option::Start{"["},
        indicators::option::Fill{"="}, indicators::option::Lead{">"},
        indicators::option::Remainder{" "}, indicators::option::End{"]"},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true},
        indicators::option::MaxProgress{total},
        indicators::option::ForegroundColor{indicators::Color::cyan});
        
    // Force initial display
    bar_ptr->set_option(indicators::option::PrefixText{
        "Sample 0/" + std::to_string(total) + " "});
    bar_ptr->set_progress(0);
    std::cerr.flush();
  }

  if (bar_ptr) {
    // Update every 5% or every 10 samples, whichever is smaller
    std::size_t update_interval = std::min<std::size_t>(total / 20, 10);
    update_interval = std::max<std::size_t>(update_interval, 1);
    
    bool should_update = (i == 0) || (i >= total) || 
                        ((i - last_update) >= update_interval) ||
                        ((i % update_interval) == 0);
    
    if (should_update) {
      // Update the prefix to show the current sample count.
      bar_ptr->set_option(indicators::option::PrefixText{
          "Sample " + std::to_string(i) + "/" + std::to_string(total) + " "});

      bar_ptr->set_progress(i);
      std::cerr.flush();
      last_update = i;
    }

    // When the loop is finished, mark as complete and clean up.
    if (i >= total) {
      // Clear cyan bar first, then show completion message
      bar_ptr->mark_as_completed();
      std::cerr.flush();
      
      // Reset the pointer to ensure a new bar is created for any subsequent fits
      bar_ptr.reset();
    }
  }
}

/**
 * @brief Main entry point for chi-square fitting procedure
 * @param chisq_ref Reference to ChiSquare object (e.g., SpectrumFit)
 * @param csm_info Minimizer configuration information
 * @param chisq_dof Output chi-square per degree of freedom
 * @param fitqual Output fit quality (p-value)
 * @param bestfit_params Output best-fit parameter estimates
 * @param param_covariance Output parameter covariance matrix
 * @param out_sampling_file Output file for sampling results
 * @param xmlout XML handler for output
 * @param kobs Observable handler for data management
 *
 * Automatically detects MPI environment and chooses appropriate execution mode:
 * - Single process: Serial execution
 * - Multiple processes: MPI parallel execution
 */
void doChiSquareFitting(ChiSquare& chisq_ref,
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual,
                        vector<MCEstimate>& bestfit_params,
                        RealSymmetricMatrix& param_covariance,
                        const std::string& out_sampling_file,
                        XMLHandler& xmlout, KBObsHandler* kobs) {

  // Check if we're running with multiple MPI processes
  int current_size, current_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &current_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

  // If we have multiple processes, use traditional MPI mode
  if (current_size > 1) {
    if (current_rank == 0) {
      std::cerr << "Running ChiSquareFitting with " << current_size
                << "MPI processes" << std::endl;
    }

    // Call the traditional MPI version that uses MPI_COMM_WORLD directly
    doChiSquareFittingMPI(chisq_ref, csm_info, chisq_dof, fitqual,
                          bestfit_params, param_covariance, out_sampling_file,
                          xmlout, kobs, MPI_COMM_WORLD);
    return;
  }
  std::cerr << "Running ChiSquareFitting in serial mode" << std::endl;
  doChiSquareFittingSerial(chisq_ref, csm_info, chisq_dof, fitqual,
                           bestfit_params, param_covariance, out_sampling_file,
                           xmlout, kobs);
}

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
 * Performs the complete fitting workflow in serial mode:
 * 1. Full sample fit to get initial best-fit parameters
 * 2. Bootstrap/jackknife resampling fits for uncertainty estimation
 * 3. Result aggregation and output generation
 *
 * The progress is tracked with a visual progress bar.
 */
void doChiSquareFittingSerial(ChiSquare& chisq_ref,
                              const ChiSquareMinimizerInfo& csm_info,
                              double& chisq_dof, double& fitqual,
                              vector<MCEstimate>& bestfit_params,
                              RealSymmetricMatrix& param_covariance,
                              const std::string& out_sampling_file,
                              XMLHandler& xmlout, KBObsHandler* kobs) {
  uint nparams = chisq_ref.getNumberOfFitParameters();
  uint nresiduals = chisq_ref.getNumberOfResiduals();
  uint nsamplings = chisq_ref.getNumberOfResamplings();
  // TODO: switch back to 2
  if (nresiduals <= (nparams + 1))
    throw(std::invalid_argument("Too few residuals, fit cannot proceed"));
  uint dof = nresiduals - nparams;
  MCEnsembleInfo mcindep(kobs->getNumberOfResamplings());

  // assign initial guess for fit parameters
  uint sampindex = 0; // full sample
  chisq_ref.setResamplingIndex(sampindex);

  vector<double> start_params;
  start_params.resize(nparams);
  bool only_update_prior_initial_guesses = false;
  chisq_ref.guessInitialFitParamValues(start_params,
                                       only_update_prior_initial_guesses);
  vector<MCObsInfo> fitparaminfos;
  fitparaminfos.resize(nparams);
  chisq_ref.getFitParamMCObsInfo(fitparaminfos);

  vector<KBObsInfo> kbfitparaminfos;
  set<KBObsInfo> kbfitparaminfoset;
  for (vector<MCObsInfo>::const_iterator it = fitparaminfos.begin();
       it != fitparaminfos.end(); ++it) {
    kbfitparaminfos.push_back(KBObsInfo(mcindep, *it));
    kbfitparaminfoset.insert(KBObsInfo(mcindep, *it));
  }

  stringstream logger;
  logger << "Degrees of freedom  = " << dof << endl;

  // set up the minimizer
  ChiSquareMinimizer CSM(chisq_ref, csm_info);

  vector<double> params_fullsample;
  XMLHandler xmlz;
  double chisq;
  RealSymmetricMatrix pcov;

  // minimize using the full sampling

  std::cerr << "Starting minimization with full sample" << std::endl;

  bool flag =
      CSM.findMinimum(start_params, chisq, params_fullsample, pcov, xmlz);

  if (xmlz.good())
    xmlout.put_child(xmlz);
  if (!flag) {
    throw(std::invalid_argument("Fitting with full sample failed"));
  }
  chisq_dof = chisq / double(dof);
  std::cerr << "Full sample fit completed with chi-square = " << chisq
            << " and chi-square per dof = " << chisq_dof << std::endl;
  logger << "Full sample chisq/dof = " << chisq_dof << endl;
  for (uint p = 0; p < nparams; ++p) {
    logger << "params_fullsample[" << p << "] = " << params_fullsample[p]
           << endl;
  }

  for (uint p = 0; p < nparams; ++p)
    kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_fullsample[p]);

  vector<double> params_sample;

  //   loop over the re-samplings
  list<uint> failed;
  char origverbose = CSM.getVerbosity();
  CSM.setVerbosity('L'); // quiet inner fits

  const std::size_t N = nsamplings; // total samples
  only_update_prior_initial_guesses = true;

  std::cerr << "Starting minimization with" << N << " resamplings" << std::endl;
  for (sampindex = 1; sampindex <= N; ++sampindex) {
    vector<double> start(params_fullsample);
    double chisq_samp;
    chisq_ref.setResamplingIndex(sampindex);
    chisq_ref.guessInitialFitParamValues(start,
                                         only_update_prior_initial_guesses);

    bool flag = CSM.findMinimum(start, chisq_samp, params_sample);

    // progress bar (handled internally by Indicators)
    show_progress(sampindex, N);

    // detailed per-sample log
    logger << "Resamplings index = " << sampindex << " chisq = " << chisq_samp
           << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = " << params_sample[p] << '\n';

    if (flag) {
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_sample[p]);
    } else {
      logger << "Above fit failed!\n";
      failed.push_back(sampindex);
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_fullsample[p]);
    }
  }

  std::cerr << '\n';             // finish the progress line
  CSM.setVerbosity(origverbose); // restore caller’s verbosity

  xmlformat("ResamplingsMinimizationsLog", logger.str(), xmlz);
  if (xmlz.good())
    xmlout.put_child(xmlz);

  XMLHandler xmlso;
  kobs->writeSamplingValuesToFile(kbfitparaminfoset, out_sampling_file, xmlso,
                                  true);
  xmlout.put_child(xmlso);

  bestfit_params.resize(nparams);
  XMLHandler xmlres("BestFitResult");
  xmlres.put_child("NumberOfResiduals", make_string(nresiduals));
  xmlres.put_child("NumberOfFitParameters", make_string(nparams));
  xmlres.put_child("DegreesOfFreedom", make_string(dof));
  xmlres.put_child("ChiSquarePerDof", make_string(chisq_dof));
  fitqual = getChiSquareFitQuality(dof, chisq);
  xmlres.put_child("FitQuality", make_string(fitqual));
  for (uint p = 0; p < nparams; ++p) {
    XMLHandler xmlp("FitParameter" + make_string(p));
    XMLHandler xmlpi;
    fitparaminfos[p].output(xmlpi);
    xmlp.put_child(xmlpi);
    // if param has a name, output it
    string obs_name = fitparaminfos[p].getObsName();
    auto& registry = ParameterNameRegistry::getInstance();
    string param_name = registry.getParameterNameFromMCObsName(obs_name);
    if (!param_name.empty()) {
      xmlp.put_child("ParameterName", param_name);
    }
    bestfit_params[p] = kobs->getEstimate(kbfitparaminfos[p]);
    XMLHandler xmlfp;
    bestfit_params[p].output(xmlfp);
    xmlp.put_child(xmlfp);
    xmlres.put_child(xmlp);
  }

  XMLHandler xmlcov("FitParameterCovariances");
  for (uint p = 0; p < nparams; ++p) {
    for (uint pp = p; pp < nparams; ++pp) {
      double cov = kobs->getCovariance(kbfitparaminfos[p], kbfitparaminfos[pp]);
      xmlcov.put_child("Cov_" + make_string(p) + "_" + make_string(pp),
                       make_string(cov));
    }
  }
  xmlres.put_child(xmlcov);
  xmlout.put_child(xmlres);
}

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
 * Performs chi-square fitting with MPI parallelization:
 * 1. Rank 0 performs full sample fit
 * 2. Results are broadcast to all ranks
 * 3. Resampling fits are distributed across ranks in round-robin fashion
 * 4. Results are gathered back to rank 0
 * 5. Rank 0 aggregates results and generates output
 *
 * Features:
 * - Load-balanced work distribution
 * - Efficient result aggregation with MPI_Gatherv
 * - Progress tracking on rank 0
 * - Comprehensive error handling and logging
 */
void doChiSquareFittingMPI(
    ChiSquare& chisq_ref, const ChiSquareMinimizerInfo& csm_info,
    double& chisq_dof, double& fitqual, vector<MCEstimate>& bestfit_params,
    RealSymmetricMatrix& param_covariance, const std::string& out_sampling_file,
    XMLHandler& xmlout, KBObsHandler* kobs, MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // check if using slurm/srun
  bool is_slurm = false;
  if (rank == 0) {
    char* slurm_job_id = getenv("SLURM_JOB_ID");
    if (slurm_job_id != nullptr) {
      is_slurm = true;
      std::cerr << "Running in SLURM environment with job ID: " << slurm_job_id
                    << std::endl;
    }
  }

  uint nparams = chisq_ref.getNumberOfFitParameters();
  uint nresiduals = chisq_ref.getNumberOfResiduals();
  uint nsamplings = chisq_ref.getNumberOfResamplings();

  if (nresiduals <= (nparams + 1))
    throw(std::invalid_argument("Too few residuals, fit cannot proceed"));
  uint dof = nresiduals - nparams;

  MCEnsembleInfo mcindep(kobs->getNumberOfResamplings());

  vector<double> start_params;
  start_params.resize(nparams);
  vector<MCObsInfo> fitparaminfos;
  fitparaminfos.resize(nparams);
  chisq_ref.getFitParamMCObsInfo(fitparaminfos);

  vector<KBObsInfo> kbfitparaminfos;
  set<KBObsInfo> kbfitparaminfoset;
  for (vector<MCObsInfo>::const_iterator it = fitparaminfos.begin();
       it != fitparaminfos.end(); ++it) {
    kbfitparaminfos.push_back(KBObsInfo(mcindep, *it));
    kbfitparaminfoset.insert(KBObsInfo(mcindep, *it));
       }

  stringstream logger;
  vector<double> params_fullsample;
  double chisq = 0.0;
  RealSymmetricMatrix pcov;

  // Full sample fit (only rank 0)
  if (rank == 0) {
    try {
      logger << "Degrees of freedom  = " << dof << endl;

      // Initial guess for fit parameters
      uint sampindex = 0; // full sample
      chisq_ref.setResamplingIndex(sampindex);
      bool only_update_prior_initial_guesses = false;
      chisq_ref.guessInitialFitParamValues(start_params,
                                           only_update_prior_initial_guesses);

      // Set up the minimizer
      ChiSquareMinimizer CSM(chisq_ref, csm_info);

      XMLHandler xmlz;
      std::cerr << "Starting minimization with full sample" << std::endl;

      bool flag =
          CSM.findMinimum(start_params, chisq, params_fullsample, pcov, xmlz);

      if (xmlz.good())
        xmlout.put_child(xmlz);
      if (!flag) {
        throw(std::runtime_error("Full sample fit failed"));
      }

      chisq_dof = chisq / double(dof);
      std::cerr << "Full sample fit completed with chi-square = " << chisq
                << " and chi-square per dof = " << chisq_dof << std::endl;
      logger << "Full sample chisq/dof = " << chisq_dof << endl;
      for (uint p = 0; p < nparams; ++p) {
        logger << "params_fullsample[" << p << "] = " << params_fullsample[p]
               << endl;
      }

      // Store full sample results
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_fullsample[p]);

    } catch (const std::exception& e) {
      std::cerr << "Full sample fit failed: " << e.what() << std::endl;
      throw;
    }
  }

  // Ensure all ranks have properly sized vectors before broadcast
  if (rank != 0) {
    params_fullsample.resize(nparams);
  }

  // Broadcast the full-sample fit parameters to all ranks (workers only need
  // the parameters)
  MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, comm);

  // Ensure all ranks have consistent data before starting resamplings
  MPI_Barrier(comm);

  if (rank == 0) {
    std::cerr << "Starting parallel minimization with " << nsamplings << " resamplings across "
              << size << " MPI ranks" << std::endl;
  }

  // Each rank processes its assigned resamplings
  ChiSquareMinimizer CSM(chisq_ref, csm_info);
  char origverbose = CSM.getVerbosity();
  CSM.setVerbosity('L'); // quiet inner fits

  vector<double> params_sample;
  bool only_update_prior_initial_guesses = true;

  // Track failed fits for logging (local list per rank)
  list<uint> failed_local;

  // Distribute samples among ranks (round-robin)
  vector<uint> my_samples;
  for (uint sampindex = 1; sampindex <= nsamplings; ++sampindex) {
    if ((sampindex - 1) % size == static_cast<uint>(rank)) {
      my_samples.push_back(sampindex);
    }
  }

  // Buffer to store this rank's fit results for communication
  std::vector<double> local_results(my_samples.size() * nparams, 0.0);
  std::size_t local_pos = 0;

  auto t0 = std::chrono::steady_clock::now();

  // Initialize counters for the rank-0 progress bar
  std::size_t completed_local = 0; // only for estimating global progress
  // progress bar handled inside show_progress; no extra bookkeeping needed
  const std::size_t N = nsamplings; // total samples
  only_update_prior_initial_guesses = true;

  // Process assigned resamplings
  for (uint sampindex : my_samples) {

    vector<double> start(params_fullsample);
    double chisq_samp;

    // Set up ChiSquare for this sample
    chisq_ref.setResamplingIndex(sampindex);
    chisq_ref.guessInitialFitParamValues(start,
                                         only_update_prior_initial_guesses);

    // Perform minimization (start vector will be modified and carry forward)
    bool flag = CSM.findMinimum(start, chisq_samp, params_sample);

    // Log results (all ranks log)
    logger << "Rank " << rank << " Resamplings index = " << sampindex
           << " chisq = " << chisq_samp << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = " << params_sample[p] << '\n';

    // Store results in local buffer for later communication
    if (!flag) {
      logger << "Above fit failed!\n";
      failed_local.push_back(sampindex);
      params_sample = params_fullsample; // fall back to full-sample values
    }
    for (uint p = 0; p < nparams; ++p)
      local_results[local_pos++] = params_sample[p];

    // Show progress (rank 0 only). We extrapolate the work done by rank 0
    if (rank == 0) {
      ++completed_local; // one more sample finished on rank 0
      std::size_t global_est =
          std::min<std::size_t>(completed_local * size, nsamplings);
      show_progress(global_est, nsamplings, is_slurm);
    }
  }

  CSM.setVerbosity(origverbose);

  // Critical: Ensure ALL ranks have completed their samples
  // before rank 0 tries to access the results
  MPI_Barrier(comm);

  if (rank == 0) {
    std::cerr << "\nAll ranks completed resampling fits. Aggregating results..."
              << std::endl;
    std::cerr.flush();
  }

  // ------------------------------------------------------------
  // Gather fit results from all ranks to rank 0
  // ------------------------------------------------------------
  int sendcount = static_cast<int>(local_results.size());
  std::vector<int> recvcounts;
  std::vector<int> displs;
  std::vector<double> all_results;
  MPI_Gather(&sendcount, 1, MPI_INT,
             rank == 0 ? (recvcounts.resize(size), recvcounts.data()) : nullptr,
             1, MPI_INT, 0, comm);

  if (rank == 0) {
    displs.resize(size, 0);
    int total = 0;
    for (int r = 0; r < size; ++r) {
      displs[r] = total;
      total += recvcounts[r];
    }
    all_results.resize(total);
  }

  MPI_Gatherv(local_results.data(), sendcount, MPI_DOUBLE,
              rank == 0 ? all_results.data() : nullptr,
              rank == 0 ? recvcounts.data() : nullptr,
              rank == 0 ? displs.data() : nullptr, MPI_DOUBLE, 0, comm);

  // Gather failed sample indices
  std::vector<uint> failed_vec(failed_local.begin(), failed_local.end());
  int failed_count = static_cast<int>(failed_vec.size());
  std::vector<int> recvcounts_fail;
  std::vector<int> displs_fail;
  std::vector<uint> all_failed;
  MPI_Gather(&failed_count, 1, MPI_INT,
             rank == 0 ? (recvcounts_fail.resize(size), recvcounts_fail.data())
                       : nullptr,
             1, MPI_INT, 0, comm);

  if (rank == 0) {
    displs_fail.resize(size, 0);
    int total_fail = 0;
    for (int r = 0; r < size; ++r) {
      displs_fail[r] = total_fail;
      total_fail += recvcounts_fail[r];
    }
    all_failed.resize(total_fail);
  }

  MPI_Gatherv(failed_vec.data(), failed_count, MPI_UNSIGNED,
              rank == 0 ? all_failed.data() : nullptr,
              rank == 0 ? recvcounts_fail.data() : nullptr,
              rank == 0 ? displs_fail.data() : nullptr, MPI_UNSIGNED, 0, comm);

  // Rank 0 inserts gathered results into its KBObsHandler
  if (rank == 0) {
    auto compute_samples = [&](int r) {
      std::vector<uint> s;
      for (uint samp = 1; samp <= nsamplings; ++samp)
        if ((samp - 1) % size == static_cast<uint>(r))
          s.push_back(samp);
      return s;
    };

    for (int r = 0; r < size; ++r) {
      auto samp_r = compute_samples(r);
      for (std::size_t idx = 0; idx < samp_r.size(); ++idx) {
        uint sampindex = samp_r[idx];
        std::size_t base = displs[r] + idx * nparams;
        for (uint p = 0; p < nparams; ++p) {
          double val = all_results[base + p];
          kobs->putSamplingValue(kbfitparaminfos[p], sampindex, val);
        }
      }
    }

    for (uint fi : all_failed) {
      failed_local.push_back(fi);
    }
  }

  MPI_Barrier(comm);

  // Collect logger strings from all ranks for combined output
  if (rank == 0) {
    std::cerr << "Aggregating logs from all ranks..."
              << std::endl;
    // Start with rank 0's logger
    std::ostringstream combined_logger;
    combined_logger << logger.str();

    // Collect logger strings from other ranks
    for (int r = 1; r < size; ++r) {
      int log_size;
      MPI_Recv(&log_size, 1, MPI_INT, r, 100, comm, MPI_STATUS_IGNORE);
      std::vector<char> log_buffer(log_size + 1);
      MPI_Recv(log_buffer.data(), log_size, MPI_CHAR, r, 101, comm,
               MPI_STATUS_IGNORE);
      log_buffer[log_size] = '\0';

      // Append this rank's log to combined output
      combined_logger << log_buffer.data();
    }

    XMLHandler xmlz;
    xmlformat("ResamplingsMinimizationsLog", combined_logger.str(), xmlz);
    if (xmlz.good())
      xmlout.put_child(xmlz);
  } else {
    // Send logger string to rank 0
    std::string my_log = logger.str();
    int log_size = static_cast<int>(my_log.size());
    MPI_Send(&log_size, 1, MPI_INT, 0, 100, comm);
    MPI_Send(my_log.c_str(), log_size, MPI_CHAR, 0, 101, comm);
  }

  // Additional barrier to ensure kobs writes are fully synchronized
  MPI_Barrier(comm);

  // Only rank 0 generates final output (all sampling data should now be in
  // kobs)
  if (rank == 0) {
    // Generate remaining output (same structure as serial version)

    XMLHandler xmlso;
    kobs->writeSamplingValuesToFile(kbfitparaminfoset, out_sampling_file, xmlso,
                                    true);
    xmlout.put_child(xmlso);

    bestfit_params.resize(nparams);
    XMLHandler xmlres("BestFitResult");
    xmlres.put_child("NumberOfResiduals", make_string(nresiduals));
    xmlres.put_child("NumberOfFitParameters", make_string(nparams));
    xmlres.put_child("DegreesOfFreedom", make_string(dof));
    xmlres.put_child("ChiSquarePerDof", make_string(chisq_dof));
    fitqual = getChiSquareFitQuality(dof, chisq);
    xmlres.put_child("FitQuality", make_string(fitqual));

    for (uint p = 0; p < nparams; ++p) {
      XMLHandler xmlp("FitParameter" + make_string(p));
      XMLHandler xmlpi;
      fitparaminfos[p].output(xmlpi);
      xmlp.put_child(xmlpi);
      string obs_name = fitparaminfos[p].getObsName();
      auto& registry = ParameterNameRegistry::getInstance();
      string param_name = registry.getParameterNameFromMCObsName(obs_name);
      if (!param_name.empty()) {
        xmlp.put_child("ParameterName", param_name);
      }
      bestfit_params[p] = kobs->getEstimate(kbfitparaminfos[p]);
      XMLHandler xmlfp;
      bestfit_params[p].output(xmlfp);
      xmlp.put_child(xmlfp);
      xmlres.put_child(xmlp);
    }

    XMLHandler xmlcov("FitParameterCovariances");
    for (uint p = 0; p < nparams; ++p) {
      for (uint pp = p; pp < nparams; ++pp) {
        double cov =
            kobs->getCovariance(kbfitparaminfos[p], kbfitparaminfos[pp]);
        xmlcov.put_child("Cov_" + make_string(p) + "_" + make_string(pp),
                         make_string(cov));
      }
    }
    xmlres.put_child(xmlcov);
    xmlout.put_child(xmlres);

    if (!failed_local.empty()) {
      XMLHandler xmlf("FailedResamplings");
      for (list<uint>::const_iterator it = failed_local.begin();
           it != failed_local.end(); ++it) {
        xmlf.put_child("FailedIndex", make_string(*it));
           }
      xmlout.put_child(xmlf);
    }
  }

  // Final synchronization
  MPI_Barrier(comm);
  if (rank == 0) {
    std::cerr << "Chi-square fitting completed successfully." << std::endl;
  }
}

/**
 * @name Gamma Function Utilities
 * @brief Incomplete gamma function implementations for chi-square statistics
 *
 * These functions provide implementations of the incomplete gamma function
 * and related utilities needed for computing chi-square fit quality (p-values).
 *
 * Functions provided:
 * - gammln(xx): Natural logarithm of gamma function
 * - Qgamma(s, x): Upper incomplete gamma function ratio
 * - Pgamma(s, x): Lower incomplete gamma function ratio
 * - getChiSquareFitQuality(dof, chisq): Chi-square p-value calculation
 *
 * The implementation uses series expansion for small x and continued fraction
 * for large x to ensure numerical stability across the full domain.
 *
 * @{
 */

/**
 * @brief Natural logarithm of the gamma function
 * @param xx Input value (must be > 0)
 * @return ln(Gamma(xx))
 *
 * Uses Lanczos approximation for numerical stability.
 */
double gammln(double xx) {
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146,     -86.50532032941677,
                          24.01409824083091,     -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

void gcf(double& gammcf, double a, double x, double& gln) {
  int i;
  double an, b, c, d, del, h;
  const double eps = 3.0e-12;
  const int itmax = 100;
  const double fpmin = 1.0e-30;

  gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0 / fpmin;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= itmax; i++) {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < fpmin)
      d = fpmin;
    c = b + an / c;
    if (fabs(c) < fpmin)
      c = fpmin;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < eps)
      break;
  }
  if (i > itmax) {
    throw(std::invalid_argument("a too large, itmax too small in gcf"));
  }
  gammcf = exp(-x + a * log(x) - gln) * h;
}

void gser(double& gamser, double a, double x, double& gln) {
  int n;
  double sum, del, ap;
  const int itmax = 100;
  const double eps = 3.0e-12;

  gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)
      throw(std::invalid_argument("x less than 0 in routine gser"));
    gamser = 0.0;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= itmax; n++) {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs(del) < fabs(sum) * eps) {
        gamser = sum * exp(-x + a * log(x) - gln);
        return;
      }
    }
    throw(
        std::invalid_argument("a too large, itmax too small in routine gser"));
  }
}

double Qgamma(double a, double x) {
  double gamser, gammcf, gln;
  if (x < 0.0 || a <= 0.0) {
    throw(std::invalid_argument("Invalid arguments in routine gammq"));
  }
  if (x < (a + 1.0)) {
    gser(gamser, a, x, gln);
    return 1.0 - gamser;
  } else {
    gcf(gammcf, a, x, gln);
    return gammcf;
  }
}

double Pgamma(double a, double x) {
  double gamser, gammcf, gln;
  if (x < 0.0 || a <= 0.0) {
    throw(std::invalid_argument("Invalid arguments in routine gammp"));
  }
  if (x < (a + 1.0)) {
    gser(gamser, a, x, gln);
    return gamser;
  } else {
    gcf(gammcf, a, x, gln);
    return 1.0 - gammcf;
  }
}

/**
 * @brief Calculate chi-square fit quality (p-value)
 * @param dof Degrees of freedom
 * @param chisquare Chi-square value
 * @return Fit quality (p-value)
 *
 * Returns the chi-square quality of fit, given by:
 *     Q = Qgamma(dof/2, chisquare/2)
 *
 * This is the probability that a chi-square random variable with
 * 'dof' degrees of freedom would exceed the observed value.
 *
 * Interpretation:
 * - Values near 1: Very good fit
 * - Values near 0: Very poor fit
 * - Values < 0.05: Typically considered poor fit
 */
double getChiSquareFitQuality(unsigned int dof, double chisquare) {
  return Qgamma(0.5 * double(dof), 0.5 * chisquare);
}

/** @} */ // End of Gamma Function Utilities

// ****************************************************************************
