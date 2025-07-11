#include "chisq_fit.h"
#include <mpi.h>
#include <sstream>
#include <algorithm>
#include <thread>
using namespace std;

// *************************************************************************

// Serial version (original implementation)
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

  vector<double> start(params_fullsample);
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
// MPI version using traditional MPI launch (srun/mpirun)
// *************************************************************************

void doChiSquareFitting(ChiSquare& chisq_ref,
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual,
                        vector<MCEstimate>& bestfit_params,
                        RealSymmetricMatrix& param_covariance,
                        const std::string& out_sampling_file,
                        XMLHandler& xmlout, KBObsHandler* kobs) {
  
  // Check if we're running with multiple MPI processes
  int current_size;
  MPI_Comm_size(MPI_COMM_WORLD, &current_size);
  
  // If we have multiple processes, use traditional MPI mode
  if (current_size > 1) {
    std::cerr << "Using traditional MPI mode with " << current_size << " processes for chi-square fitting" << std::endl;
    
    // Call the traditional MPI version that uses MPI_COMM_WORLD directly
    doChiSquareFittingMPI_Traditional(chisq_ref, csm_info, chisq_dof, fitqual,
                                     bestfit_params, param_covariance, out_sampling_file,
                                     xmlout, kobs, MPI_COMM_WORLD);
    return;
  }
  
  // Single process mode - use serial execution
  std::cerr << "Using serial execution" << std::endl;
  doChiSquareFittingSerial(chisq_ref, csm_info, chisq_dof, fitqual,
                          bestfit_params, param_covariance, out_sampling_file,
                          xmlout, kobs);
}



// *************************************************************************
// Traditional MPI version using existing processes (original approach)
// *************************************************************************

void doChiSquareFittingMPI_Traditional(ChiSquare& chisq_ref,
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
  
  // Debug: Show rank assignment
  if (rank == 0) {
    std::cerr << "Distributing work across " << size << " MPI ranks" << std::endl;
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
  int fit_success = 1; // 1 = success, 0 = failure
  
  // Only rank 0 does the full sample fit
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
        fit_success = 0;
      } else {
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
      }
    } catch (const std::exception& e) {
      fit_success = 0;
      std::cerr << "Full sample fit failed: " << e.what() << std::endl;
    }
  }
  
  // Broadcast fit success status to all ranks
  MPI_Bcast(&fit_success, 1, MPI_INT, 0, comm);
  if (!fit_success) {
    throw(std::invalid_argument("Fitting with full sample failed"));
  }
  
  // Broadcast full sample results to all ranks
  if (rank == 0) {
    start_params = params_fullsample;
  }
  
  // Ensure all ranks have properly sized arrays before broadcast
  start_params.resize(nparams);
  params_fullsample.resize(nparams);
  
  // Broadcast the results
  MPI_Bcast(&chisq_dof, 1, MPI_DOUBLE, 0, comm);
  MPI_Bcast(start_params.data(), nparams, MPI_DOUBLE, 0, comm);
  MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, comm);
  
  // Distribute resamplings among ranks
  vector<uint> my_samples;
  for (uint sampindex = 1; sampindex <= nsamplings; ++sampindex) {
    if ((sampindex - 1) % size == (uint)rank) {
      my_samples.push_back(sampindex);
    }
  }
  
  // Debug: Show sample distribution (only rank 0 for summary)
  if (rank == 0) {
    std::cerr << "Each rank will process approximately " << (nsamplings + size - 1) / size << " samples" << std::endl;
  }
  
  // Each rank processes its assigned resamplings
  ChiSquareMinimizer CSM(chisq_ref, csm_info);
  char origverbose = CSM.getVerbosity();
  CSM.setVerbosity('L'); // quiet inner fits
  
  vector<double> params_sample;
  bool only_update_prior_initial_guesses = true;
  
  // Store results locally instead of writing to shared kobs
  vector<pair<uint, vector<double>>> my_successful_results;
  vector<pair<uint, vector<double>>> my_failed_results;
  list<uint> my_failed;
  
  if (rank == 0) {
    std::cerr << "Starting minimization with resamplings across " << size << " MPI ranks" << std::endl;
  }
  
  // For progress tracking (rank 0 only)
  auto t0 = std::chrono::steady_clock::now();
  int last_percent = -1; // Track progress percentage changes
  
  // Process assigned resamplings
  for (size_t i = 0; i < my_samples.size(); ++i) {
    uint sampindex = my_samples[i];
    double chisq_samp;
    
    chisq_ref.setResamplingIndex(sampindex);
    vector<double> start(start_params);
    chisq_ref.guessInitialFitParamValues(start, only_update_prior_initial_guesses);
    
    bool flag = CSM.findMinimum(start, chisq_samp, params_sample);
    
    // Show progress (rank 0 only) - similar to serial version
    if (rank == 0) {
      // Estimate total progress based on rank 0's work
      uint estimated_total = (i + 1) * size;
      if (estimated_total > nsamplings) estimated_total = nsamplings;
      
      // Show progress when percentage changes (like the serial version)
      int percent = static_cast<int>(100.0 * estimated_total / nsamplings + 0.5);
      if (percent != last_percent || estimated_total == nsamplings) {
        show_progress(estimated_total, nsamplings, t0, std::cerr);
        std::cerr.flush(); // Ensure output is visible
        last_percent = percent;
      }
    }
    
    // Detailed per-sample log
    logger << "Rank " << rank << " Resamplings index = " << sampindex
           << " chisq = " << chisq_samp << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = " << params_sample[p] << '\n';
    
    // Store results locally instead of directly writing to kobs
    if (flag) {
      my_successful_results.push_back(make_pair(sampindex, params_sample));
    } else {
      logger << "Above fit failed!\n";
      my_failed.push_back(sampindex);
      my_failed_results.push_back(make_pair(sampindex, params_fullsample));
    }
  }
  
  CSM.setVerbosity(origverbose);
  
  // Wait for all ranks to complete
  MPI_Barrier(comm);
  
  if (rank == 0) {
    std::cerr << '\n'; // finish the progress line
    std::cerr << "Gathering results from all ranks..." << std::endl;
    
    // Write rank 0's results to kobs first
    for (const auto& result : my_successful_results) {
      uint sampindex = result.first;
      const vector<double>& params = result.second;
      for (uint p = 0; p < nparams; ++p) {
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params[p]);
      }
    }
    for (const auto& result : my_failed_results) {
      uint sampindex = result.first;
      const vector<double>& params = result.second;
      for (uint p = 0; p < nparams; ++p) {
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params[p]);
      }
    }
    
    // Gather failed samples for output
    list<uint> all_failed(my_failed);
    
    // Gather results from other ranks
    for (int r = 1; r < size; ++r) {
      // Receive successful results
      int num_successful;
      MPI_Recv(&num_successful, 1, MPI_INT, r, 0, comm, MPI_STATUS_IGNORE);
      for (int i = 0; i < num_successful; ++i) {
        int sampindex;
        vector<double> params(nparams);
        MPI_Recv(&sampindex, 1, MPI_INT, r, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(params.data(), nparams, MPI_DOUBLE, r, 2, comm, MPI_STATUS_IGNORE);
        
        for (uint p = 0; p < nparams; ++p) {
          kobs->putSamplingValue(kbfitparaminfos[p], static_cast<uint>(sampindex), params[p]);
        }
      }
      
      // Receive failed results
      int num_failed;
      MPI_Recv(&num_failed, 1, MPI_INT, r, 3, comm, MPI_STATUS_IGNORE);
      for (int i = 0; i < num_failed; ++i) {
        int sampindex;
        vector<double> params(nparams);
        MPI_Recv(&sampindex, 1, MPI_INT, r, 4, comm, MPI_STATUS_IGNORE);
        MPI_Recv(params.data(), nparams, MPI_DOUBLE, r, 5, comm, MPI_STATUS_IGNORE);
        
        for (uint p = 0; p < nparams; ++p) {
          kobs->putSamplingValue(kbfitparaminfos[p], static_cast<uint>(sampindex), params[p]);
        }
        all_failed.push_back(static_cast<uint>(sampindex));
      }
    }
    
    // Generate output
    XMLHandler xmlz;
    xmlformat("ResamplingsMinimizationsLog", logger.str(), xmlz);
    if (xmlz.good())
      xmlout.put_child(xmlz);
    
    XMLHandler xmlso;
    kobs->writeSamplingValuesToFile(kbfitparaminfoset, out_sampling_file, xmlso, true);
    xmlout.put_child(xmlso);
    
    bestfit_params.resize(nparams);
    XMLHandler xmlres("BestFitResult");
    xmlres.put_child("NumberOfResiduals", make_string(nresiduals));
    xmlres.put_child("NumberOfFitParameters", make_string(nparams));
    xmlres.put_child("DegreesOfFreedom", make_string(dof));
    xmlres.put_child("ChiSquarePerDof", make_string(chisq_dof));
    double fitqual_local = getChiSquareFitQuality(dof, chisq);
    fitqual = fitqual_local;
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
    
    if (all_failed.size() > 0) {
      XMLHandler xmlf("FailedResamplings");
      for (list<uint>::const_iterator it = all_failed.begin(); it != all_failed.end(); ++it) {
        xmlf.put_child("FailedIndex", make_string(*it));
      }
      xmlout.put_child(xmlf);
    }
    
  } else {
    // Non-root ranks send their results
    
    // Send successful results
    int num_successful = static_cast<int>(my_successful_results.size());
    MPI_Send(&num_successful, 1, MPI_INT, 0, 0, comm);
    for (const auto& result : my_successful_results) {
      int sampindex = static_cast<int>(result.first);
      MPI_Send(&sampindex, 1, MPI_INT, 0, 1, comm);
      MPI_Send(result.second.data(), nparams, MPI_DOUBLE, 0, 2, comm);
    }
    
    // Send failed results
    int num_failed = static_cast<int>(my_failed_results.size());
    MPI_Send(&num_failed, 1, MPI_INT, 0, 3, comm);
    for (const auto& result : my_failed_results) {
      int sampindex = static_cast<int>(result.first);
      MPI_Send(&sampindex, 1, MPI_INT, 0, 4, comm);
      MPI_Send(result.second.data(), nparams, MPI_DOUBLE, 0, 5, comm);
    }
  }
  
  // Final barrier to ensure all ranks complete together
  MPI_Barrier(comm);
  
  // Broadcast final results to all ranks
  MPI_Bcast(&fitqual, 1, MPI_DOUBLE, 0, comm);
  
  // All ranks need to have the final bestfit_params size
  if (rank != 0) {
    bestfit_params.resize(nparams);
  }
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
