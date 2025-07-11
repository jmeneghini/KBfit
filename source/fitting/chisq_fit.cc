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
// MPI version that spawns processes dynamically
// *************************************************************************

void doChiSquareFitting(ChiSquare& chisq_ref,
                        const ChiSquareMinimizerInfo& csm_info,
                        double& chisq_dof, double& fitqual,
                        vector<MCEstimate>& bestfit_params,
                        RealSymmetricMatrix& param_covariance,
                        const std::string& out_sampling_file,
                        XMLHandler& xmlout, KBObsHandler* kobs,
                        int num_mpi_processes) {
  
  // If only 1 process requested, use serial version
  if (num_mpi_processes <= 1) {
    doChiSquareFittingSerial(chisq_ref, csm_info, chisq_dof, fitqual,
                            bestfit_params, param_covariance, out_sampling_file,
                            xmlout, kobs);
    return;
  }
  
  std::cerr << "Starting dynamic MPI chi-square fitting with " << num_mpi_processes << " processes" << std::endl;
  
  uint nparams = chisq_ref.getNumberOfFitParameters();
  uint nresiduals = chisq_ref.getNumberOfResiduals();
  uint nsamplings = chisq_ref.getNumberOfResamplings();
  
  if (nresiduals <= (nparams + 1))
    throw(std::invalid_argument("Too few residuals, fit cannot proceed"));
  uint dof = nresiduals - nparams;
  
  MCEnsembleInfo mcindep(kobs->getNumberOfResamplings());

  // Initial guess for fit parameters
  vector<double> start_params;
  start_params.resize(nparams);
  uint sampindex = 0; // full sample
  chisq_ref.setResamplingIndex(sampindex);
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

  // Set up the minimizer and do full sample fit
  ChiSquareMinimizer CSM(chisq_ref, csm_info);
  vector<double> params_fullsample;
  XMLHandler xmlz;
  double chisq;
  RealSymmetricMatrix pcov;

  std::cerr << "Starting minimization with full sample" << std::endl;
  bool flag = CSM.findMinimum(start_params, chisq, params_fullsample, pcov, xmlz);

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
    logger << "params_fullsample[" << p << "] = " << params_fullsample[p] << endl;
  }

  // Store full sample results
  for (uint p = 0; p < nparams; ++p)
    kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_fullsample[p]);

  // Now spawn MPI processes for parallel resampling work
  MPI_Comm intercomm;
  int spawn_result;
  
  // Get the current executable name and arguments to spawn worker processes
  char* executable = nullptr;
  char** argv = nullptr;
  
  // For simplicity, we'll use a self-spawning approach where we pass a special flag
  // to indicate this is a worker process
  char worker_flag[] = "--mpi-worker";
  char* worker_argv[] = {worker_flag, nullptr};
  
  std::cerr << "Spawning " << num_mpi_processes << " MPI worker processes..." << std::endl;
  
  // Spawn worker processes
  spawn_result = MPI_Comm_spawn("KBfit", worker_argv, num_mpi_processes,
                                MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, 
                                MPI_ERRCODES_IGNORE);
  
  if (spawn_result != MPI_SUCCESS) {
    std::cerr << "Warning: MPI_Comm_spawn failed, falling back to serial execution" << std::endl;
    // Fall back to serial processing
    doChiSquareFittingSerial(chisq_ref, csm_info, chisq_dof, fitqual,
                            bestfit_params, param_covariance, out_sampling_file,
                            xmlout, kobs);
    return;
  }

  // Create merged communicator
  MPI_Comm merged_comm;
  MPI_Intercomm_merge(intercomm, 0, &merged_comm);  // 0 = low group (parent)
  
  // Now distribute the resampling work using the merged communicator
  int merged_rank, merged_size;
  MPI_Comm_rank(merged_comm, &merged_rank);
  MPI_Comm_size(merged_comm, &merged_size);
  
  std::cerr << "Merged communicator created with " << merged_size << " processes" << std::endl;
  
  // Broadcast necessary data to all processes
  // Send basic parameters
  uint params_data[4] = {nparams, nresiduals, nsamplings, dof};
  MPI_Bcast(params_data, 4, MPI_UNSIGNED, 0, merged_comm);
  
  // Send full sample results
  MPI_Bcast(&chisq_dof, 1, MPI_DOUBLE, 0, merged_comm);
  MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, merged_comm);
  
  // Send minimizer info (serialize to string)
  XMLHandler csm_xml;
  csm_info.output(csm_xml);
  string csm_str = csm_xml.output();
  int csm_str_len = csm_str.length();
  MPI_Bcast(&csm_str_len, 1, MPI_INT, 0, merged_comm);
  MPI_Bcast(const_cast<char*>(csm_str.c_str()), csm_str_len, MPI_CHAR, 0, merged_comm);
  
  // TODO: We would need to serialize and send the ChiSquare object and KBObsHandler
  // This is complex and depends on the specific implementation
  // For now, we'll fall back to serial execution
  
  std::cerr << "Warning: Full MPI implementation requires serialization of complex objects." << std::endl;
  std::cerr << "Falling back to serial execution for now." << std::endl;
  
  // Clean up MPI communicators
  MPI_Comm_free(&merged_comm);
  MPI_Comm_free(&intercomm);
  
  // Fall back to serial processing of resamplings
  vector<double> start(params_fullsample);
  vector<double> params_sample;
  list<uint> failed;
  char origverbose = CSM.getVerbosity();
  CSM.setVerbosity('L');
  only_update_prior_initial_guesses = true;

  auto t0 = std::chrono::steady_clock::now();
  std::cerr << "Starting minimization with resamplings (serial fallback)" << std::endl;
  
  for (sampindex = 1; sampindex <= nsamplings; ++sampindex) {
    double chisq_samp;
    chisq_ref.setResamplingIndex(sampindex);
    vector<double> start_samp(start);
    chisq_ref.guessInitialFitParamValues(start_samp, only_update_prior_initial_guesses);
    
    bool flag = CSM.findMinimum(start_samp, chisq_samp, params_sample);
    
    // Show progress
    if (sampindex % max(1u, nsamplings / 20) == 0) {
      show_progress(sampindex, nsamplings, t0, std::cerr);
    }
    
    logger << "Resamplings index = " << sampindex << " chisq = " << chisq_samp << '\n';
    for (uint p = 0; p < nparams; ++p)
      logger << "params_sample[" << p << "] = " << params_sample[p] << '\n';
    
    if (flag) {
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_sample[p]);
    } else {
      logger << "Above fit failed!\n";
      failed.push_back(sampindex);
      for (uint p = 0; p < nparams; ++p)
        kobs->putSamplingValue(kbfitparaminfos[p], sampindex, params_fullsample[p]);
    }
  }
  
  CSM.setVerbosity(origverbose);
  std::cerr << '\n'; // finish the progress line
  
  // Generate output (same as serial version)
  xmlformat("ResamplingsMinimizationsLog", logger.str(), xmlz);
  if (xmlz.good())
    xmlout.put_child(xmlz);
  
  if (!out_sampling_file.empty()) {
    XMLHandler xmlso;
    kobs->writeSamplingValuesToFile(kbfitparaminfoset, out_sampling_file, xmlso, true);
    xmlout.put_child(xmlso);
  }

  if (failed.size() > 0) {
    XMLHandler xmlf("FailedResamplings");
    for (list<uint>::const_iterator it = failed.begin(); it != failed.end(); ++it) {
      xmlf.put_child("FailedIndex", make_string(*it));
    }
    xmlout.put_child(xmlf);
  }

  // Results
  bestfit_params.resize(nparams);
  param_covariance.resize(nparams);
  
  XMLHandler xmlres("FitResults");
  xmlres.put_child("NumberOfFitParameters", make_string(nparams));
  xmlres.put_child("NumberOfResiduals", make_string(nresiduals));
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
}

// *************************************************************************
// MPI worker function for spawned processes
// *************************************************************************

void doChiSquareFittingMPI(ChiSquare& chisq_ref,
                           const ChiSquareMinimizerInfo& csm_info,
                           double& chisq_dof, double& fitqual,
                           vector<MCEstimate>& bestfit_params,
                           RealSymmetricMatrix& param_covariance,
                           const std::string& out_sampling_file,
                           XMLHandler& xmlout, KBObsHandler* kobs,
                           int num_processes) {
  
  // This function will be called by spawned MPI processes
  // It receives work from the master process and processes assigned resamplings
  
  MPI_Comm parent_comm;
  MPI_Comm_get_parent(&parent_comm);
  
  if (parent_comm == MPI_COMM_NULL) {
    std::cerr << "Error: MPI worker called but no parent communicator found" << std::endl;
    return;
  }
  
  // Merge with parent communicator
  MPI_Comm merged_comm;
  MPI_Intercomm_merge(parent_comm, 1, &merged_comm);  // 1 = high group (child)
  
  int rank, size;
  MPI_Comm_rank(merged_comm, &rank);
  MPI_Comm_size(merged_comm, &size);
  
  std::cerr << "MPI worker process " << rank << " of " << size << " started" << std::endl;
  
  // Receive basic parameters
  uint params_data[4];
  MPI_Bcast(params_data, 4, MPI_UNSIGNED, 0, merged_comm);
  uint nparams = params_data[0];
  uint nresiduals = params_data[1]; 
  uint nsamplings = params_data[2];
  uint dof = params_data[3];
  
  // Receive full sample results
  double chisq_dof_recv;
  MPI_Bcast(&chisq_dof_recv, 1, MPI_DOUBLE, 0, merged_comm);
  
  vector<double> params_fullsample(nparams);
  MPI_Bcast(params_fullsample.data(), nparams, MPI_DOUBLE, 0, merged_comm);
  
  // Receive minimizer info
  int csm_str_len;
  MPI_Bcast(&csm_str_len, 1, MPI_INT, 0, merged_comm);
  string csm_str(csm_str_len, '\0');
  MPI_Bcast(const_cast<char*>(csm_str.c_str()), csm_str_len, MPI_CHAR, 0, merged_comm);
  
  std::cerr << "MPI worker " << rank << " received parameters and ready to process resamplings" << std::endl;
  
  // TODO: Receive and reconstruct ChiSquare and KBObsHandler objects
  // This requires complex serialization
  
  // For now, just participate in the communicator and then exit
  MPI_Comm_free(&merged_comm);
  
  std::cerr << "MPI worker " << rank << " finished (serialization not implemented)" << std::endl;
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
