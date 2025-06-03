#include "chisq_spectrum.h"
using namespace std;

SpectrumFit::SpectrumFit(XMLHandler& xmlin,

}



SpectrumFit::~SpectrumFit() {
  clear();
}

void SpectrumFit::clear() {
  for (uint k = 0; k < BQ.size(); ++k)
    delete BQ[k];
  delete Kmat;
  delete Kinv;
  BQ.clear();
}

// Need to add non-Kmatrix fit parameters below here
// (L, mass, etc.), just use mean values.
void SpectrumFit::guessInitialFitParamValues(
    vector<double>& fitparams) const {
  if (Kmat != 0)
    fitparams = Kmat->getParameterValues();
  else
    fitparams = Kinv->getParameterValues();
}

void SpectrumFit::getFitParamMCObsInfo(
    vector<MCObsInfo>& fitinfos) const {
  const vector<KFitParamInfo>* fpptr = 0;
  if (Kmat != 0)
    fpptr = &(Kmat->getFitParameterInfos());
  else
    fpptr = &(Kinv->getFitParameterInfos());
  fitinfos.resize(fpptr->size());
  for (uint k = 0; k < fitinfos.size(); ++k) {
    fitinfos[k] = MCObsInfo((*fpptr)[k].getMCObsName());
  }
}

// might want to modify this at some point once I
// better understand its purpose
void SpectrumFit::do_output(XMLHandler& xmlout) const {
  xmlout.set_root("DeterminantResidualFit");
  XMLHandler xmlK;
  if (Kmat != 0)
    Kmat->output(xmlK);
  else
    Kinv->output(xmlK);
  xmlout.put_child(xmlK);
}

void SpectrumFit::evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) {

}


// By introducing the priors for our MC observables,
// cov(r_i, r_j) simplifies to cov(R_i, R_j), where
// R is the observable
void SpectrumFit::initializeInvCov() {
  bool bootstrap = KBOH->isBootstrapMode();

  // get params
  vector<MCObsInfo> param_infos;
  getFitParamMCObsInfo(param_infos);
  vector<vector<double>> obs_samples(param_infos.size());

  RealSymmetricMatrix cov(nresiduals, 0.0);
  for (uint k = 0; k < nresiduals; ++k)
    for (uint j = 0; j <= k; ++j) {
      if ((j == k) || (ensemble_id[k] == ensemble_id[j]))
        cov(k, j) = (bootstrap) ? KBOH->boot_covariance(res[k], res[j])
                                : KBOH->jack_covariance(res[k], res[j]);
    }
}

