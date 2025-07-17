#include "minimizer.h"
#include <iostream>
using namespace std;

// ***************************************************************************

ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(XMLHandler& xmlin) {
  XMLHandler xmlr(xmlin, "MinimizerInfo");
  m_method = defaultmethod;
  string reply;
  if (xmlreadifchild(xmlr, "Method", reply)) {
    if (reply == "NL2Sno")
      m_method = 'N';
#ifndef NO_MINUIT
    else if (reply == "Minuit2Migrad")
      m_method = 'M';
    else if (reply == "Minuit2Simplex")
      m_method = 'S';
#endif
    else
      throw(std::invalid_argument(
          "Invalid <Method> tag in ChiSquareMinimizerInfo"));
  }
  m_param_reltol = 1e-6;
  double tol;
  if (xmlreadifchild(xmlr, "ParameterRelTol", tol))
    setParameterRelativeTolerance(tol);
  m_chisq_reltol = 1e-4;
  if (xmlreadifchild(xmlr, "ChiSquareRelTol", tol))
    setChiSquareRelativeTolerance(tol);
  m_max_its = 1024;
  int its;
  if (xmlreadifchild(xmlr, "MaximumIterations", its))
    setMaximumIterations(its);
  m_verbosity = 'L';
  if (xmlreadifchild(xmlr, "Verbosity", reply)) {
    if (reply == "Low")
      m_verbosity = 'L';
    else if (reply == "Medium")
      m_verbosity = 'M';
    else if (reply == "High")
      m_verbosity = 'H';
    else
      throw(std::invalid_argument(
          "Invalid <Verbosity> tag in ChiSquareMinimizerInfo"));
  }
}

ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(char method, double param_reltol,
                                               double chisq_reltol, int max_its,
                                               char verbosity) {
  setMethod(method);
  setParameterRelativeTolerance(param_reltol);
  setChiSquareRelativeTolerance(chisq_reltol);
  setMaximumIterations(max_its);
  setVerbosity(verbosity);
}

ChiSquareMinimizerInfo::ChiSquareMinimizerInfo(
    const ChiSquareMinimizerInfo& info)
    : m_method(info.m_method), m_param_reltol(info.m_param_reltol),
      m_chisq_reltol(info.m_chisq_reltol), m_max_its(info.m_max_its),
      m_verbosity(info.m_verbosity) {}

ChiSquareMinimizerInfo&
ChiSquareMinimizerInfo::operator=(const ChiSquareMinimizerInfo& info) {
  m_method = info.m_method;
  m_param_reltol = info.m_param_reltol;
  m_chisq_reltol = info.m_chisq_reltol;
  m_max_its = info.m_max_its;
  m_verbosity = info.m_verbosity;
  return *this;
}

void ChiSquareMinimizerInfo::setMethod(char method) {
  if ((method == 'N')
#ifndef NO_MINUIT
      || (method == 'M') || (method == 'S')
#endif
  ) {
    m_method = method;
    return;
  }
  m_method = defaultmethod;
  throw(std::invalid_argument(
      "Invalid method character in ChiSquareMinimizerInfo::setMethod"));
}

void ChiSquareMinimizerInfo::setMinuit2Simplex() {
#ifndef NO_MINUIT
  m_method = 'S';
#else
  throw(std::invalid_argument(
      "Minuit2 library not available in ChiSquareMinimizerInfo::setMethod"));
#endif
}

void ChiSquareMinimizerInfo::setMinuit2Migrad() {
#ifndef NO_MINUIT
  m_method = 'M';
#else
  throw(std::invalid_argument(
      "Minuit2 library not available in ChiSquareMinimizerInfo::setMethod"));
#endif
}

void ChiSquareMinimizerInfo::setParameterRelativeTolerance(double rtol) {
  if (rtol > 0.0) {
    m_param_reltol = rtol;
    return;
  }
  throw(std::invalid_argument(
      "Invalid input in "
      "ChiSquareMinimizerInfo::setParameterRelativeTolerance"));
}

void ChiSquareMinimizerInfo::setChiSquareRelativeTolerance(double rtol) {
  if (rtol > 0.0) {
    m_chisq_reltol = rtol;
    return;
  }
  throw(std::invalid_argument(
      "Invalid input in "
      "ChiSquareMinimizerInfo::setChiSquareRelativeTolerance"));
}

void ChiSquareMinimizerInfo::setMaximumIterations(unsigned int maxit) {
  if (maxit > 10) {
    m_max_its = maxit;
    return;
  }
  throw(std::invalid_argument(
      "Invalid input in ChiSquareMinimizerInfo::setMaximumIterations"));
}

void ChiSquareMinimizerInfo::setVerbosity(char verbosity) {
  if ((verbosity == 'L') || (verbosity == 'M') || (verbosity == 'H')) {
    m_verbosity = verbosity;
    return;
  }
  throw(std::invalid_argument(
      "Invalid verbosity character in ChiSquareMinimizerInfo::setVerbosity"));
}

void ChiSquareMinimizerInfo::output(XMLHandler& xmlout) const {
  xmlout.set_root("MinimizerInfo");
  if (m_method == 'S')
    xmlout.put_child("Method", "Minuit2Simplex");
  else if (m_method == 'M')
    xmlout.put_child("Method", "Minuit2Migrad");
  else if (m_method == 'N')
    xmlout.put_child("Method", "NL2Sno");
  xmlout.put_child("ParameterRelTol", make_string(m_param_reltol));
  xmlout.put_child("ChiSquareRelTol", make_string(m_chisq_reltol));
  xmlout.put_child("MaximumIterations", make_string(m_max_its));
  if (m_verbosity == 'L')
    xmlout.put_child("Verbosity", "Low");
  else if (m_verbosity == 'M')
    xmlout.put_child("Verbosity", "Medium");
  else if (m_verbosity == 'H')
    xmlout.put_child("Verbosity", "High");
}

string ChiSquareMinimizerInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string ChiSquareMinimizerInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

// ***************************************************************************

ChiSquareMinimizer::ChiSquareMinimizer(ChiSquare& in_chisq)
    : m_chisq(&in_chisq), m_nl2sno(0)
#ifndef NO_MINUIT
      ,
      m_minuit2(0)
#endif
{
  alloc_method();
}

ChiSquareMinimizer::ChiSquareMinimizer(ChiSquare& in_chisq,
                                       const ChiSquareMinimizerInfo& info)
    : m_chisq(&in_chisq), m_info(info), m_nl2sno(0)
#ifndef NO_MINUIT
      ,
      m_minuit2(0)
#endif
{
  alloc_method();
}

ChiSquareMinimizer::~ChiSquareMinimizer() { dealloc_method(); }

void ChiSquareMinimizer::dealloc_method() {
  delete m_nl2sno;
  m_nl2sno = 0;
#ifndef NO_MINUIT
  delete m_minuit2;
  m_minuit2 = 0;
#endif
}

void ChiSquareMinimizer::alloc_method() {
  if (m_info.m_method == 'N')
    m_nl2sno = new NL2SnoMinimizer(*m_chisq);
#ifndef NO_MINUIT
  else if ((m_info.m_method == 'S') || (m_info.m_method == 'M'))
    m_minuit2 = new Minuit2ChiSquare(*m_chisq);
#endif
}

void ChiSquareMinimizer::reset(const ChiSquareMinimizerInfo& info) {
  if (m_info.m_method == info.m_method)
    return;
  dealloc_method();
  m_info = info;
  alloc_method();
}

bool ChiSquareMinimizer::find_minimum(const vector<double>& starting_params,
                                      double& chisq_min,
                                      vector<double>& params_at_minimum,
                                      RealSymmetricMatrix& params_covariance,
                                      XMLHandler& xmlout, char verbosity) {
  if (m_info.m_method == 'N')
    return find_minimum_nl2sno(starting_params, chisq_min, params_at_minimum,
                               params_covariance, xmlout, verbosity);
#ifndef NO_MINUIT
  else if ((m_info.m_method == 'S') || (m_info.m_method == 'M'))
    return find_minimum_minuit2(starting_params, m_info.m_method, chisq_min,
                                params_at_minimum, params_covariance, xmlout,
                                verbosity);
#endif
  return false;
}

bool ChiSquareMinimizer::findMinimum(const vector<double>& starting_params,
                                     double& chisq_min,
                                     vector<double>& params_at_minimum,
                                     RealSymmetricMatrix& params_covariance,
                                     XMLHandler& xmlout) {
  return find_minimum(starting_params, chisq_min, params_at_minimum,
                      params_covariance, xmlout, m_info.m_verbosity);
}

bool ChiSquareMinimizer::findMinimum(double& chisq_min,
                                     vector<double>& params_at_minimum,
                                     RealSymmetricMatrix& params_covariance,
                                     XMLHandler& xmlout) {
  vector<double> starting_params(m_chisq->getNumberOfFitParameters());
  m_chisq->guessInitialFitParamValues(starting_params, false);
  return find_minimum(starting_params, chisq_min, params_at_minimum,
                      params_covariance, xmlout, m_info.m_verbosity);
}

bool ChiSquareMinimizer::findMinimum(const vector<double>& starting_params,
                                     double& chisq_min,
                                     vector<double>& params_at_minimum,
                                     RealSymmetricMatrix& params_covariance) {
  XMLHandler xmlout;
  return find_minimum(starting_params, chisq_min, params_at_minimum,
                      params_covariance, xmlout, m_info.m_verbosity);
}

bool ChiSquareMinimizer::findMinimum(double& chisq_min,
                                     vector<double>& params_at_minimum,
                                     RealSymmetricMatrix& params_covariance) {
  vector<double> starting_params(m_chisq->getNumberOfFitParameters());
  m_chisq->guessInitialFitParamValues(starting_params, false);
  XMLHandler xmlout;
  return find_minimum(starting_params, chisq_min, params_at_minimum,
                      params_covariance, xmlout, m_info.m_verbosity);
}

bool ChiSquareMinimizer::findMinimum(const vector<double>& starting_params,
                                     double& chisq_min,
                                     vector<double>& params_at_minimum) {
  XMLHandler xmlout;
  RealSymmetricMatrix params_covariance;
  return find_minimum(starting_params, chisq_min, params_at_minimum,
                      params_covariance, xmlout, m_info.m_verbosity);
}

#ifndef NO_MINUIT

bool ChiSquareMinimizer::find_minimum_minuit2(
    const vector<double>& starting_params, char method, double& chisq_min,
    vector<double>& params_at_minimum, RealSymmetricMatrix& params_covariance,
    XMLHandler& xmlout, char verbosity) {
  xmlout.clear();
  uint nparam = m_chisq->getNumberOfFitParameters();
  if (starting_params.size() != nparam)
    throw(std::invalid_argument("Invalid starting parameters"));

  /* --- OPTIONAL: turn live logging on/off from XML -------------- */
  const bool want_trace = (verbosity != 'L'); // Medium or High

  if (want_trace) {
    // 0 = silent, 1 = low, 2 = medium, 3 = high (lots of lines)
    ROOT::Minuit2::MnPrint::SetGlobalLevel(3);
  } else {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(0);
  }

  std::vector<double> unc(nparam);
  for (uint p = 0; p < nparam; ++p)
    unc[p] = 0.01 * starting_params[p]; // set up initial uncertainties

  unsigned int strategylevel = 2; // 0 = low, 1 = med, 2 = high quality
                                  // lower level means faster, higher means
                                  // more reliable minimization

  ROOT::Minuit2::MnUserParameters upar;
  for (uint p = 0; p < nparam; ++p) {
    upar.Add(string("param") + make_string(p), starting_params[p], unc[p]);
  }
  ROOT::Minuit2::MnUserParameterState param(upar);
  ROOT::Minuit2::MnStrategy strategy(strategylevel);

  ROOT::Minuit2::FunctionMinimum* csmin;
  if (method == 'M') {
    //    ROOT::Minuit2::MnMigrad M(*m_minuit2, starting_params, unc,
    //    strategylevel);
    ROOT::Minuit2::MnMigrad M(*m_minuit2, param, strategy);
    csmin = new ROOT::Minuit2::FunctionMinimum(
        M(m_info.m_max_its, m_info.m_chisq_reltol));
  } else {
    //    ROOT::Minuit2::MnSimplex M(*m_minuit2, starting_params, unc,
    //    strategylevel);
    ROOT::Minuit2::MnSimplex M(*m_minuit2, param, strategy);
    csmin = new ROOT::Minuit2::FunctionMinimum(
        M(m_info.m_max_its, m_info.m_chisq_reltol));
  }

  if (verbosity != 'L') {
    ostringstream outlog;
    outlog << "Minuit2 Minimization Result:\n" << (*csmin);
    xmlformat("Minuit2Log", outlog.str(), xmlout);
  }

  bool flag = csmin->IsValid();
  chisq_min = csmin->Fval();
  params_at_minimum.resize(nparam);
  for (uint p = 0; p < nparam; ++p)
    params_at_minimum[p] = csmin->UserParameters().Value(p);
  if (flag) {
    if (method == 'M') {
      params_covariance.resize(nparam);
      for (uint p = 0; p < nparam; ++p)
        for (uint q = p; q < nparam; ++q)
          params_covariance(p, q) = csmin->UserCovariance()(p, q);
    } else {
      params_covariance.clear();
    }
  }
  delete csmin;

  return flag;
}

#endif

bool ChiSquareMinimizer::find_minimum_nl2sno(
    const vector<double>& starting_params, double& chisq_min,
    vector<double>& params_at_minimum, RealSymmetricMatrix& param_cov,
    XMLHandler& xmlout, char verbosity) {
  xmlout.clear();
  ostringstream outlog;
  m_nl2sno->setInitialFitParamValues(starting_params);
  int flag = m_nl2sno->chisq_fit(m_info.m_param_reltol, m_info.m_chisq_reltol,
                                 m_info.m_verbosity, m_info.m_max_its,
                                 chisq_min, outlog);

  if (verbosity != 'L')
    xmlformat("NL2SnoLog", outlog.str(), xmlout);

  params_at_minimum = m_nl2sno->getFitParams();
  if (flag <= 2) {
    param_cov = m_nl2sno->getParamCovariance();
    return true;
  }

  param_cov.clear();
  return false;
}

// ********************************************************************

#ifndef NO_MINUIT

void Minuit2ChiSquare::guessInitialFitParamValues(vector<double>& params) {
  m_chisq->guessInitialFitParamValues(params, false);
}

double Minuit2ChiSquare::operator()(const vector<double>& params) const {
  return m_chisq->evalChiSquare(params);
}

double Minuit2ChiSquare::evalChiSquare(const vector<double>& params) const {
  return m_chisq->evalChiSquare(params);
}

#endif

// ********************************************************************

namespace MinPack {
int nl2sno(ChiSquare* M, int n, int p, vector<double>& params,
           RealSymmetricMatrix& param_cov, double& chisq,
           vector<double>& residuals, vector<int>& viv, vector<double>& vv,
           ostringstream& outlog);
int dfault(int* iv, double* v);
}; // namespace MinPack

int NL2SnoMinimizer::chisq_fit(double paramreltol, double chisqreltol,
                               char verbosity, int max_its, double& chisq,
                               ostringstream& outlog) {
  int p = m_chisq->getNumberOfFitParameters();
  int n = m_chisq->getNumberOfResiduals();

  //  set up nl2sno parameters

  MinPack::dfault(&iv[0], &v[0]);
  iv[14] = 1;           // print covariance matrix
  iv[16] = 4 * max_its; // max func evals
  iv[17] = max_its;     // max iterations
  iv[18] =
      (verbosity == 'L') ? 0 : ((verbosity == 'H') ? 1 : 10); // output flag
  v[31] = chisqreltol; // relative func tolerance
  v[32] = paramreltol; // relative solution tolerance

  //  perform minimization

  int info = MinPack::nl2sno(m_chisq, n, p, m_fitparams, m_param_cov, chisq,
                             m_residuals, iv, v, outlog);

  //  return code

  if (iv[0] == 5)
    info = 0; // both x and f convergence
  else if (iv[0] == 3)
    info = 1; // x convergence only
  else if (iv[0] == 4)
    info = 2; // f convergence only
  else if ((iv[0] == 9) || (iv[0] == 10))
    info = 3; // max its exceeded
  else if ((iv[0] == 7) || (iv[0] == 8))
    info = 4; // false or singular convergence
  else
    info = 5; // bad input or miscellaneous

  return info;
}

void NL2SnoMinimizer::guessInitialFitParamValues() {
  m_chisq->guessInitialFitParamValues(m_fitparams, false);
}

void NL2SnoMinimizer::setInitialFitParamValues(
    const std::vector<double>& start_params) {
  if (start_params.size() != m_nparams)
    throw(std::invalid_argument("Invalid starting parameters"));
  for (uint k = 0; k < m_nparams; ++k)
    m_fitparams[k] = start_params[k];
}

// ************************************************************
