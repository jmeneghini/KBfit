#include "chisq_fit.h"
#include <mpi.h>
#include <sstream>
#include <algorithm>
#include <thread>
#include <cstdio>
using namespace std;

// *************************************************************************

// Runs the chi-square fitting procedure, either in serial or MPI mode.

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
      std::cerr << "Running ChiSquareFitting with " << current_size << "MPI processes" << std::endl;
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

// Serial version
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
  chisq_ref.guessInitialFitParamValues(start_params, only_update_prior_initial_guesses);
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
  std::cerr << "Full sample fit completed with chi-square = "
            << chisq << " and chi-square per dof = " << chisq_dof << std::endl;
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
  CSM.setVerbosity('L');                             // quiet inner fits

  auto   t0              = std::chrono::steady_clock::now();
  int    last_percent    = -1;                       // nothing printed yet
  const  std::size_t N   = nsamplings;               // total samples
  only_update_prior_initial_guesses = true;

  std::cerr << "Starting minimization with resamplings" << std::endl;
  for (sampindex = 1; sampindex <= N; ++sampindex)
  {
    vector<double> start(params_fullsample);
    double chisq_samp;
    chisq_ref.setResamplingIndex(sampindex);
    chisq_ref.guessInitialFitParamValues(start, only_update_prior_initial_guesses);

    bool flag = CSM.findMinimum(start, chisq_samp, params_sample);

    // ── progress bar ------------------------------------------------------
    int percent = static_cast<int>(100.0 * sampindex / N + 0.5); // round
    if (percent != last_percent || sampindex == N) {
      show_progress(sampindex, N, t0);          // updates the bar
      last_percent = percent;
    }

    // ── detailed per-sample log------------------------------
    logger << "Resamplings index = " << sampindex
           << " chisq = "            << chisq_samp << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = "
             << params_sample[p] << '\n';

    if (flag) {
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_sample[p]);
    } else {
      logger << "Above fit failed!\n";
      failed.push_back(sampindex);
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_fullsample[p]);
    }
  }

  std::cerr << '\n';               // finish the progress line
  CSM.setVerbosity(origverbose);   // restore caller’s verbosity

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
    string param_name =
        registry.getParameterNameFromMCObsName(obs_name);
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

// *************************************************************************
// MPI Parallel version with coordinated ChiSquare access
// *************************************************************************

void doChiSquareFittingMPI(ChiSquare& chisq_ref,
                                     const ChiSquareMinimizerInfo& csm_info,
                                     double& chisq_dof, double& fitqual,
                                     vector<MCEstimate>& bestfit_params,
                                     RealSymmetricMatrix& param_covariance,
                                     const std::string& out_sampling_file,
                                     XMLHandler& xmlout, KBObsHandler* kobs,
                                     MPI_Comm comm) {
  
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
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
      chisq_ref.guessInitialFitParamValues(start_params, only_update_prior_initial_guesses);
      
      // Set up the minimizer
      ChiSquareMinimizer CSM(chisq_ref, csm_info);
      
      XMLHandler xmlz;
      std::cerr << "Starting minimization with full sample" << std::endl;
      
      bool flag = CSM.findMinimum(start_params, chisq, params_fullsample, pcov, xmlz);
      
      if (xmlz.good())
        xmlout.put_child(xmlz);
      if (!flag) {
        throw(std::runtime_error("Full sample fit failed"));
      }
      
      chisq_dof = chisq / double(dof);
      std::cerr << "Full sample fit completed with chi-square = "
                << chisq << " and chi-square per dof = " << chisq_dof << std::endl;
      logger << "Full sample chisq/dof = " << chisq_dof << endl;
      for (uint p = 0; p < nparams; ++p) {
        logger << "params_fullsample[" << p << "] = " << params_fullsample[p] << endl;
      }
      
      // Store full sample results
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_fullsample[p]);
        
    } catch (const std::exception& e) {
      std::cerr << "Full sample fit failed: " << e.what() << std::endl;
      throw;
    }
  }

  // Ensure all ranks have properly sized vectors before broadcast
  if (rank != 0) {
    params_fullsample.resize(nparams);
  }
  
  // Broadcast the full-sample fit parameters to all ranks (workers only need the parameters)
  MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, comm);
  
  // Ensure all ranks have consistent data before starting resampling
  MPI_Barrier(comm);
  
  if (rank == 0) {
    std::cerr << "Starting parallel minimization with resamplings across " << size << " MPI ranks" << std::endl;
  }
  
  // Each rank processes its assigned resamplings
  ChiSquareMinimizer CSM(chisq_ref, csm_info);
  char origverbose = CSM.getVerbosity();
  CSM.setVerbosity('L'); // quiet inner fits
  
  vector<double> params_sample;
  bool only_update_prior_initial_guesses = true;
  
  // Track failed fits for logging
  list<uint> failed;
  

  
  // Distribute samples among ranks (round-robin)
  vector<uint> my_samples;
  for (uint sampindex = 1; sampindex <= nsamplings; ++sampindex) {
    if ((sampindex - 1) % size == static_cast<uint>(rank)) {
      my_samples.push_back(sampindex);
    }
  }
  
  auto t0 = std::chrono::steady_clock::now();
  
  // Initialize counters for the rank-0 progress bar
  std::size_t completed_local = 0;
  int         last_percent    = -1;

  // Process assigned resamplings
  for (uint sampindex : my_samples) {
    
    vector<double> start(params_fullsample);
    double chisq_samp;
    
    // Set up ChiSquare for this sample
    chisq_ref.setResamplingIndex(sampindex);
    chisq_ref.guessInitialFitParamValues(start, only_update_prior_initial_guesses);
    
    // Perform minimization (start vector will be modified and carry forward)
    bool flag = CSM.findMinimum(start, chisq_samp, params_sample);
    
    // Log results (all ranks log)
    logger << "Rank " << rank << " Resamplings index = " << sampindex
           << " chisq = " << chisq_samp << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = " << params_sample[p] << '\n';
    
    // Store results locally (don't write to kobs yet)
    if (flag) {
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_sample[p]);
    } else {
      logger << "Above fit failed!\n";
      failed.push_back(sampindex);
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex,
                               params_fullsample[p]);
    }
    
    // Show progress (rank 0 only). We extrapolate the work done by rank 0
    // to the whole set assuming balanced round-robin distribution. This
    // avoids inter-rank communication while giving a smooth progress bar
    // identical to the serial version.
    if (rank == 0) {
      ++completed_local; // one more sample finished on rank 0
      std::size_t global_est = std::min<std::size_t>(completed_local * size,
                                                     nsamplings);
      int percent = static_cast<int>(100.0 * global_est / nsamplings + 0.5);
      if (percent != last_percent || global_est == nsamplings) {
        show_progress(global_est, nsamplings, t0);
        last_percent = percent;
      }
    }
  }
  

  
  CSM.setVerbosity(origverbose);
  
  // Critical: Ensure ALL ranks have completed writing their samples to kobs
  // before rank 0 tries to access the results
  MPI_Barrier(comm);
  
  if (rank == 0) {
    std::cerr << "\nAll ranks completed resampling fits. Finalizing results..." << std::endl;
  }
  
  // ─────────────────────────────────────────────────────────────
  // (Removed per-rank sample count exchange – it was only used
  //   for debugging and generated unnecessary MPI traffic.)
  // ─────────────────────────────────────────────────────────────
  
  // Collect logger strings from all ranks for combined output
  if (rank == 0) {
    // Start with rank 0's logger
    std::ostringstream combined_logger;
    combined_logger << logger.str();
    
    // Collect logger strings from other ranks
    for (int r = 1; r < size; ++r) {
      int log_size;
      MPI_Recv(&log_size, 1, MPI_INT, r, 100, comm, MPI_STATUS_IGNORE);
      std::vector<char> log_buffer(log_size + 1);
      MPI_Recv(log_buffer.data(), log_size, MPI_CHAR, r, 101, comm, MPI_STATUS_IGNORE);
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
  
  // Only rank 0 generates final output (all sampling data should now be in kobs)
  if (rank == 0) {
    // Generate remaining output (same structure as serial version)
    
    XMLHandler xmlso;
    kobs->writeSamplingValuesToFile(kbfitparaminfoset, out_sampling_file, xmlso, true);
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
        double cov = kobs->getCovariance(kbfitparaminfos[p], kbfitparaminfos[pp]);
        xmlcov.put_child("Cov_" + make_string(p) + "_" + make_string(pp), make_string(cov));
      }
    }
    xmlres.put_child(xmlcov);
    xmlout.put_child(xmlres);
    
    if (failed.size() > 0) {
      XMLHandler xmlf("FailedResamplings");
      for (list<uint>::const_iterator it = failed.begin(); it != failed.end(); ++it) {
        xmlf.put_child("FailedIndex", make_string(*it));
      }
      xmlout.put_child(xmlf);
    }
  }
  
  // Final synchronization
  MPI_Barrier(comm);
}

// *************************************************
//
// The ratios of incomplete gamma function are defined below:
//
//       double Qgamma(double s, double x);    s>0, x>=0
//       double Pgamma(double s, double x);
//
// If an error occurs, an exception is thrown.
//
// *************************************************

// Returns the value of ln(Gamma(xx)) for xx>0

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

//  Returns the chi-square quality of fit, given by
//
//       Qgamma( dof/2,  chisquare/2 )

double getChiSquareFitQuality(unsigned int dof, double chisquare) {
  return Qgamma(0.5 * double(dof), 0.5 * chisquare);
}



// ****************************************************************************
