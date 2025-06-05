#include "chisq_base.h"
using namespace std;

// *************************************************************

void ChiSquare::initialize_base(uint number_fit_parameters,
                                uint number_residuals,
                                uint number_resamplings) {
  if (number_fit_parameters < 1)
    throw(std::invalid_argument("Number of fit parameters > 0 required"));
  if (number_residuals <= number_fit_parameters)
    throw(std::invalid_argument(
        "Number of residuals > number fit parameters required"));
  nresiduals = number_residuals;
  nfitparams = number_fit_parameters;
  inv_cov_cholesky.resize(nresiduals);
  residuals.resize(nresiduals);
  nresamplings = number_resamplings;
}

std::string ChiSquare::output(int indent) const {
  XMLHandler xmlout;
  do_output(xmlout);
  return xmlout.output(indent);
}

std::string ChiSquare::str() const {
  XMLHandler xmlout;
  do_output(xmlout);
  return xmlout.str();
}

void ChiSquare::setResamplingIndex(uint in_resampling_index) {
  if (in_resampling_index > nresamplings)
    throw(std::invalid_argument("Invalid resampling index in ChiSquare"));
  resampling_index = in_resampling_index;
}

double ChiSquare::evalChiSquare(const vector<double>& fitparams) {
  evalResidualsAndInvCovCholesky(fitparams);
  double chisq = 0.0;
  uint nres = residuals.size();
  for (int i = nres - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    chisq += tmp * tmp;
  }
  return chisq;
}

void ChiSquare::evalDiagonalResiduals(const vector<double>& fitparams,
                                      vector<double>& diag_residuals) {
  diag_residuals.resize(nresiduals);
  evalResidualsAndInvCovCholesky(fitparams);
  uint nres = residuals.size();
  for (int i = nres - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    diag_residuals[i] = tmp;
  }
}

void ChiSquare::evalDiagonalResiduals(const vector<double>& fitparams,
                                      double* diag_residuals) {
  evalResidualsAndInvCovCholesky(fitparams);
  for (int i = nresiduals - 1; i >= 0; --i) {
    double tmp = 0.0;
    for (int j = 0; j <= i; ++j)
      tmp += inv_cov_cholesky(i, j) * residuals[j];
    diag_residuals[i] = tmp;
  }
}

void ChiSquare::read_obs(XMLHandler& xmlin, const string& tag,
                                      bool get_name, MCObsInfo& obskey,
                                      set<MCObsInfo>& kset, string& name,
                                      const MCEnsembleInfo& mcens,
                                      map<KBObsInfo, double>& fixed_values) {
  try {
    XMLHandler xmlt(xmlin, tag);
    name.clear();
    if (get_name) {
      xmlread(xmlt, "Name", name, tag);
    }
    uint mcount = xmlt.count("MCObs") + xmlt.count("MCObservable");
    uint fcount = xmlt.count("FixedValue");
    if ((mcount + fcount) != 1)
      throw(std::invalid_argument("No MCObs/MCObservable or FixedValue"));
    if (mcount == 1) {
      obskey = MCObsInfo(xmlt);
      if ((obskey.isImaginaryPart()) || (obskey.isSimple()))
        throw(std::invalid_argument("MCObsInfo must be nonsimple and real"));
      if (obskey == MCObsInfo("KBScale"))
        throw(std::invalid_argument(
            "KBScale is reserved and cannot be an input MCObsInfo"));
      if (kset.find(obskey) != kset.end())
        throw(std::invalid_argument("duplicate MCObsInfo"));
      kset.insert(obskey);
    } else {
      double fixedvalue;
      xmlread(xmlt, "FixedValue", fixedvalue, tag);
      string kbname(tag);
      if (!name.empty())
        kbname += "_" + name;
      obskey = MCObsInfo(kbname);
      KBObsInfo kbkey(mcens, obskey);
      fixed_values.insert(make_pair(kbkey, fixedvalue));
    }
  } catch (const std::exception& xp) {
    string msg = string("For tag ") + tag;
    throw(std::invalid_argument(msg + string(": ") + xp.what()));
  }
}

//  This routine does the following:
//    -- searches only within an XML tag with name specified in "tag"
//    -- an <MCObs> or <MCObservable> must be encountered, then "obskey"
//        is assigned this tag; if "obskey" is already in "kset", an
//        exception is thrown, but if not, then "obskey" is inserted
//        into "kset".
//    -- "obskey" must be nonsimple and real.

void ChiSquare::read_obs(XMLHandler& xmlin, const string& tag,
                                      MCObsInfo& obskey, set<MCObsInfo>& kset) {
  try {
    XMLHandler xmlt(xmlin, tag);
    uint mcount = xmlt.count("MCObs") + xmlt.count("MCObservable");
    uint fcount = xmlt.count("FixedValue");
    if ((mcount != 1) || (fcount != 0))
      throw(std::invalid_argument(
          "No MCObs/MCObservable or disallowed FixedValue"));
    obskey = MCObsInfo(xmlt);
    if ((obskey.isImaginaryPart()) || (obskey.isSimple()))
      throw(std::invalid_argument("MCObsInfo must be nonsimple and real"));
    if (obskey == MCObsInfo("KBScale"))
      throw(std::invalid_argument(
          "KBScale is reserved and cannot be an input MCObsInfo"));
    if (kset.find(obskey) != kset.end())
      throw(std::invalid_argument("duplicate MCObsInfo"));
    kset.insert(obskey);
  } catch (const std::exception& xp) {
    string msg = string("For tag ") + tag;
    throw(std::invalid_argument(msg + string(": ") + xp.what()));
  }
}

// *************************************************************
