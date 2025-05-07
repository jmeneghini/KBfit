#include "kbobs_handler.h"
// #include <algorithm>
// #include <limits>

using namespace std;

// *********************************************************************************
// * *
// *   The class "KBObsHandler", a container class for Monte Carlo non-simple *
// *   observables, is defined in this file.  This is one of the most important
// *
// *   and central classes in "KBfit".  An object of this class is used to get
// the *
// *   Monte Carlo resamplings, store them in memory, and perform statistical *
// *   analysis on the data.  When calculations are carried out on the data or *
// *   fits are performed, objects of this class store them so that bootstrap
// and  *
// *   jackknife errors can be subsequently computed.  Objects of this class
// also  *
// *   carry out reading and writing of bootstrap/jackknife results to and from
// *
// *   files. *
// * *
// *********************************************************************************

KBObsHandler::KBObsHandler(const MCSamplingInfo& sampinfo) {
  m_in_handler = new MCObsGetHandler(sampinfo);
  m_logger << m_in_handler->getCurrentLog().str();
  m_in_handler->clearLog();
}

KBObsHandler::~KBObsHandler() { clear(); }

// if new mode != current, does complete clear

void KBObsHandler::setSamplingInfo(const MCSamplingInfo& sampinfo) {
  if ((m_in_handler->getSamplingInfo()) == sampinfo)
    return;
  m_in_handler->setSamplingInfo(sampinfo);
  clearSamplings();
}

SamplingMode KBObsHandler::getSamplingMode() const {
  return m_in_handler->getSamplingMode();
}

void KBObsHandler::connectSamplingFiles(const std::set<std::string>& sampfiles,
                                        bool verbose) {
  m_in_handler->connectSamplingFiles(sampfiles, verbose);
  m_logger << m_in_handler->getCurrentLog().str();
  m_in_handler->clearLog();
}

// if map "keys_to_keep" is empty:
//    -- keep all keys encountered for ALL ensembles encountered
// else if "keys_to_keep" DOES have an entry for an ensemble:
//    -- if entry is empty set, keep all keys for that ensemble
//    -- if entry is nonempty set, keep only those keys in the set
// else if "keys_to_keep" does NOT have an entry for an ensemble:
//    -- no keys are to be kept for that ensemble

void KBObsHandler::connectSamplingFiles(
    const std::set<std::string>& sampfiles,
    const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep,
    bool verbose) {
  m_in_handler->connectSamplingFiles(sampfiles, keys_to_keep, verbose);
  m_logger << m_in_handler->getCurrentLog().str();
  m_in_handler->clearLog();
}

void KBObsHandler::disconnectSamplingFiles(
    const std::set<std::string>& sampfiles) {
  m_in_handler->disconnectSamplingFiles(sampfiles);
  m_logger << m_in_handler->getCurrentLog().str();
  m_in_handler->clearLog();
}

void KBObsHandler::disconnectAllSamplingFiles() {
  m_in_handler->disconnectAllSamplingFiles();
  m_logger << m_in_handler->getCurrentLog().str();
  m_in_handler->clearLog();
}

bool KBObsHandler::isJackknifeMode() const {
  return m_in_handler->isJackknifeMode();
}

bool KBObsHandler::isBootstrapMode() const {
  return m_in_handler->isBootstrapMode();
}

std::set<MCEnsembleInfo> KBObsHandler::getEnsembleInfos() const {
  return m_in_handler->getEnsembleInfos();
}

const MCBinsInfo& KBObsHandler::getBinsInfo(const KBObsInfo& obskey) const {
  return m_in_handler->getBinsInfo(obskey.getMCEnsembleInfo());
}

const MCBinsInfo& KBObsHandler::getBinsInfo(const MCEnsembleInfo& mcens) const {
  return m_in_handler->getBinsInfo(mcens);
}

uint KBObsHandler::getNumberOfBins(const MCEnsembleInfo& mcens) const {
  return m_in_handler->getBinsInfo(mcens).getNumberOfBins();
}

const MCSamplingInfo& KBObsHandler::getSamplingInfo() const {
  return m_in_handler->getSamplingInfo();
}

unsigned int KBObsHandler::getNumberOfResamplings() const {
  return m_in_handler->getNumberOfResamplings();
}

void KBObsHandler::clearSamplings() { m_samples.clear(); }

void KBObsHandler::eraseSamplings(const KBObsInfo& obskey) {
  m_samples.erase(obskey);
}

void KBObsHandler::clearGetter() { m_in_handler->clear(); }

void KBObsHandler::clear() {
  m_samples.clear();
  m_in_handler->clear();
  clearLog();
}

/*
void KBObsHandler::getFileMap(XMLHandler& xmlout) const
{
 m_in_handler->getFileMap(xmlout);
}
*/

// *****************************************************************

bool KBObsHandler::queryFullAndSamplings(const KBObsInfo& obskey) {
  map<KBObsInfo, pair<RVector, uint>>::const_iterator dt =
      m_samples.find(obskey);
  if (dt != m_samples.end())
    if ((dt->second).second == (dt->second).first.size())
      return true;
  return m_in_handler->querySamplings(obskey);
}

//  return null pointer if not found

const std::pair<RVector, uint>*
KBObsHandler::get_full_and_sampling_values(const KBObsInfo& obskey) {
  map<KBObsInfo, pair<RVector, uint>>::const_iterator dt =
      m_samples.find(obskey);
  if (dt != m_samples.end()) { // cout << "FOUND in memory"<<endl;
    return &(dt->second);
  }
  RVector samples; // cout << "MUST SEARCH for in files"<<endl;
  if (m_in_handler->getSamplingsMaybe(obskey, samples)) {
    pair<map<KBObsInfo, pair<RVector, uint>>::iterator, bool> ret =
        m_samples.insert(make_pair(obskey, make_pair(samples, samples.size())));
    if (ret.second)
      return &((ret.first)->second);
  }
  return 0;
}

const RVector& KBObsHandler::getFullAndSamplingValues(const KBObsInfo& obskey) {
  const std::pair<RVector, uint>* resptr = get_full_and_sampling_values(obskey);
  if ((resptr == 0) || ((resptr->second) != resptr->first.size()))
    throw(std::runtime_error(
        string("KBObsHandler::getFullAndSamplingValues failed for ") +
        obskey.str()));
  return resptr->first;
}

// returns null pointer if not available

const RVector*
KBObsHandler::getFullAndSamplingValuesMaybe(const KBObsInfo& obskey) {
  const std::pair<RVector, uint>* resptr = get_full_and_sampling_values(obskey);
  if ((resptr == 0) || ((resptr->second) != resptr->first.size()))
    return 0;
  return &(resptr->first);
}

double KBObsHandler::getFullSampleValue(const KBObsInfo& obskey) {
  const std::pair<RVector, uint>* resptr = get_full_and_sampling_values(obskey);
  if ((resptr == 0) || (std::isnan(resptr->first[0])))
    throw(std::runtime_error(string("getSampling failed for ") + obskey.str() +
                             string(" for index = ") + make_string(0)));
  return resptr->first[0];
}

double KBObsHandler::getSamplingValue(const KBObsInfo& obskey, int samp_index) {
  const std::pair<RVector, uint>* resptr = get_full_and_sampling_values(obskey);
  if ((resptr == 0) || (samp_index < 0) ||
      (samp_index >= int(resptr->first.size())) ||
      (std::isnan(resptr->first[samp_index])))
    throw(std::runtime_error(string("getSampling failed for ") + obskey.str() +
                             string(" for index = ") +
                             make_string(samp_index)));
  return resptr->first[samp_index];
}

bool KBObsHandler::getSamplingValueMaybe(const KBObsInfo& obskey,
                                         int samp_index, double& result) {
  const std::pair<RVector, uint>* resptr = get_full_and_sampling_values(obskey);
  if ((resptr == 0) || (samp_index < 0) ||
      (samp_index >= int(resptr->first.size())) ||
      (std::isnan(resptr->first[samp_index])))
    return false;
  result = resptr->first[samp_index];
  return true;
}

/*
const RVector& KBObsHandler::put_samplings_in_memory(const MCObsInfo& obskey,
                      const RVector& samplings,
{
 m_samples.erase(obskey);
 pair<map<MCObsInfo,pair<RVector,uint> >::iterator,bool> ret;
 ret=m_samples.insert(make_pair(obskey,make_pair(samplings,samplings.size())));
 if (ret.second==false){
    throw(std::runtime_error("put samplings into memory failed"));}
 return ((ret.first)->second).first;
}


*/

void KBObsHandler::putSamplingValue(const KBObsInfo& obskey, int samp_index,
                                    double value, bool overwrite) {
  int sampling_max = getNumberOfResamplings();
  if (samp_index > sampling_max)
    throw(std::invalid_argument(
        "invalid index in KBObsHandler::putSamplingValue"));
  map<KBObsInfo, pair<RVector, uint>>::iterator dt = m_samples.find(obskey);
  if (dt != m_samples.end()) {
    double& entry = (dt->second.first)[samp_index];
    if (std::isnan(entry)) {
      entry = value;
      (dt->second.second)++;
      return;
    } else if (overwrite) {
      entry = value;
      return;
    } else {
      throw(std::invalid_argument(
          "cannot putCurrentSamplingValue since no overwrite"));
    }
  }
  RVector buffer(sampling_max + 1, std::numeric_limits<double>::quiet_NaN());
  buffer[samp_index] = value;
  m_samples.insert(make_pair(obskey, make_pair(buffer, 1)));
}

const RVector& KBObsHandler::putFullAndSamplings(const KBObsInfo& obskey,
                                                 const RVector& samps,
                                                 bool overwrite) {
  uint sampling_max = getNumberOfResamplings();
  if (samps.size() != (sampling_max + 1))
    throw(std::invalid_argument(
        "invalid sampling vector in KBObsHandler::putFullAndSamplings"));
  map<KBObsInfo, pair<RVector, uint>>::iterator dt = m_samples.find(obskey);
  if (dt != m_samples.end()) {
    if (overwrite) {
      dt->second.first = samps;
      return dt->second.first;
    } else {
      throw(std::invalid_argument(
          "cannot putFullAndSamplings since no overwrite"));
    }
  }
  pair<map<KBObsInfo, pair<RVector, uint>>::iterator, bool> pt =
      m_samples.insert(make_pair(obskey, make_pair(samps, samps.size())));
  if (!(pt.second))
    throw(std::invalid_argument(
        "putFullAndSamplings failed due to memory allocation problem"));
  return (pt.first)->second.first;
}

void KBObsHandler::putFixedValue(const KBObsInfo& obskey,
                                 const double& fixedvalue, uint nsamplings,
                                 bool overwrite) {
  map<KBObsInfo, pair<RVector, uint>>::iterator dt = m_samples.find(obskey);
  if (dt != m_samples.end()) {
    if (overwrite) {
      dt->second.first.resize(nsamplings);
      dt->second.first = fixedvalue;
      return;
    } else {
      throw(std::invalid_argument("cannot putFixedValue since no overwrite"));
    }
  }
  RVector samps(nsamplings + 1, fixedvalue);
  m_samples.insert(make_pair(obskey, make_pair(samps, samps.size())));
}

/*
void KBObsHandler::calc_simple_jack_samples(const RVector& bins, RVector&
samplings)
{
 uint nbins=bins.size();
 samplings.resize(nbins+1);  // 0 = full, 1..nbins are the jackknife samplings
 double dm=0.0;
 for (uint k=0;k<nbins;k++)
    dm+=bins[k];
 samplings[0]=dm/double(nbins);
 double rj=1.0/double(nbins-1);
 for (uint k=1;k<=nbins;k++)
    samplings[k]=(dm-bins[k-1])*rj;
}


void KBObsHandler::calc_simple_boot_samples(const RVector& bins, RVector&
samplings)
{
 if (Bptr==0)
    throw(std::runtime_error("No Bootstrapper object!"));
 uint nsamps=Bptr->getNumberOfResamplings();
 uint nbins=Bptr->getNumberOfObjects();
 if (nbins!=bins.size())
    throw(std::runtime_error("Mismatch in number of bins in bootstrapper"));
 samplings.resize(nsamps+1);   // 0 = full, 1..nsamps are the bootstrap
samplings double dm=0.0; for (uint k=0;k<nbins;k++) dm+=bins[k];
 samplings[0]=dm/double(nbins);
 for (uint bootindex=0;bootindex<nsamps;bootindex++){
    const Vector<uint>& indmap=Bptr->getResampling(bootindex);
    double dm=0.0;
    for (uint k=0;k<nbins;k++)
       dm+=bins[indmap[k]];
    samplings[bootindex+1]=dm/double(nbins);}
}

*/
// ***************************************************************************

MCEstimate KBObsHandler::getEstimate(const KBObsInfo& obskey) {
  if (isJackknifeMode()) {
    MCEstimate result(Jackknife);
    const RVector& sampvals = getFullAndSamplingValues(obskey);
    jack_analyze(sampvals, result);
    return result;
  } else {
    MCEstimate result(Bootstrap);
    RVector sampvals(getFullAndSamplingValues(obskey));
    boot_analyze(sampvals, result);
    return result;
  }
}

/*
MCEstimate KBObsHandler::getEstimate(const MCObsInfo& obskey, SamplingMode
inmode)
{
 MCEstimate result(inmode);
 if (inmode==Jackknife){
    const RVector& sampvals=getFullAndSamplingValues(obskey,Jackknife);
    jack_analyze(sampvals,result);
    return result;}
 else{
    RVector sampvals;
    getFullAndSamplingValues(obskey,sampvals,Bootstrap); // sampvals will be
sorted boot_analyze(sampvals,result); return result;}
}


MCEstimate KBObsHandler::getJackknifeEstimate(const MCObsInfo& obskey)
{
 MCEstimate result(Jackknife);
 const RVector& sampvals=getFullAndSamplingValues(obskey,Jackknife);
 jack_analyze(sampvals,result);
 return result;
}


MCEstimate KBObsHandler::getBootstrapEstimate(const MCObsInfo& obskey)
{
 MCEstimate result(Bootstrap);
 RVector sampvals;
 getFullAndSamplingValues(obskey,sampvals,Bootstrap); // sampvals will be sorted
 boot_analyze(sampvals,result);
 return result;
}
*/

// **********************************************************************

double KBObsHandler::getCovariance(const KBObsInfo& obskey1,
                                   const KBObsInfo& obskey2) {
  if (obskey1.getMCEnsembleInfo() != obskey2.getMCEnsembleInfo())
    return 0.0;
  const RVector& sampvals1 = getFullAndSamplingValues(obskey1);
  const RVector& sampvals2 = getFullAndSamplingValues(obskey2);
  if (sampvals1.size() != sampvals2.size())
    throw(std::runtime_error("Could not compute covariance"));
  if (isJackknifeMode())
    return jack_covariance(sampvals1, sampvals2);
  else
    return boot_covariance(sampvals1, sampvals2);
}

double KBObsHandler::getStandardError(const KBObsInfo& obskey) {
  return sqrt(getCovariance(obskey, obskey));
}

/*
double KBObsHandler::getJackKnifeError(const MCObsInfo& obskey, uint jacksize)
{
 uint nbins=getNumberOfBins();
 if ((jacksize<1)||(jacksize>(nbins/12)))
    throw(std::invalid_argument("Invalid jack size in getJackKnifeError"));
 try{
    const RVector& buffer=getBins(obskey);
    uint N=nbins - (nbins % jacksize);   // N is divisible by jacksize
    uint NJ=N/jacksize;
    double zJ=((double)(N-jacksize))/((double) N);
    double sum=buffer[0];
    for (unsigned int k=1;k<N;k++)
       sum+=buffer[k];
    double avg=sum/double(N);
    double cov=0.0;
    uint count=0;
    for (uint jack=0;jack<NJ;jack++){
       double avgJ=sum;
       for (uint i=0;i<jacksize;i++){
          avgJ-=buffer[count++];}
       avgJ/=double(N-jacksize);
       cov+=(avg-avgJ)*(avg-avgJ);}
    cov*=zJ;
    return sqrt(cov);}
 catch(const std::exception& errmsg){
    cout << "Error in KBObsHandler::getJackKnifeError: "<<errmsg.what()<<endl;
    throw;}
}
*/

//  compute jackknife error using samplings (sampvals[0] contains full sample)

void KBObsHandler::jack_analyze(const RVector& sampvals, MCEstimate& result) {
  uint n = sampvals.size() - 1;
  double avg = 0.0;
  for (uint k = 1; k <= n; k++)
    avg += sampvals[k];
  avg /= double(n); // should be the same as sampvals[0]
  double var = 0.0;
  for (uint k = 1; k <= n; k++) {
    double x = sampvals[k] - avg;
    var += x * x;
  }
  var *= (1.0 - 1.0 / double(n)); // jackknife
  result.jackassign(sampvals[0], avg, sqrt(var));
}

//  compute jackknife covariance using samplings
//  (sampvals1[0], sampvals2[0] contain full samples)

double KBObsHandler::jack_covariance(const RVector& sampvals1,
                                     const RVector& sampvals2) {
  uint n = sampvals1.size() - 1;
  double avg1 = 0.0;
  double avg2 = 0.0;
  for (uint k = 1; k <= n; k++) {
    avg1 += sampvals1[k];
    avg2 += sampvals2[k];
  }
  avg1 /= double(n);
  avg2 /= double(n);
  double jackcov = 0.0;
  for (uint k = 1; k <= n; k++) {
    jackcov += (sampvals1[k] - avg1) * (sampvals2[k] - avg2);
  }
  jackcov *= (1.0 - 1.0 / double(n)); // jackknife
  return jackcov;
}

//  compute bootstrap covariance using samplings
//  (sampvals1[0], sampvals2[0] contain full samples)

double KBObsHandler::boot_covariance(const RVector& sampvals1,
                                     const RVector& sampvals2) {
  uint n = sampvals1.size() - 1;
  double avg1 = 0.0;
  double avg2 = 0.0;
  for (uint k = 1; k <= n; k++) {
    avg1 += sampvals1[k];
    avg2 += sampvals2[k];
  }
  avg1 /= double(n);
  avg2 /= double(n);
  double bootcov = 0.0;
  for (uint k = 1; k <= n; k++) {
    bootcov += (sampvals1[k] - avg1) * (sampvals2[k] - avg2);
  }
  bootcov /= double(n - 1); // bootstrap
  return bootcov;
}

// This takes an array of "Nboot" bootstrap samples
// and returns their average (ans.boot_avg), the
// value which 84% of the samples lie above (ans.boot_low),
// the value which 84% of the samples lie below (ans.boot_upp),
// and the value which 50% of the samples lie above and
// 50% lie below (ans.boot_med).   CAUTION: this changes "sampvals"

void KBObsHandler::boot_analyze(RVector& sampvals, MCEstimate& result) {
  uint nb = sampvals.size() - 1;
  double avg = 0.0;
  for (uint k = 1; k <= nb; k++)
    avg += sampvals[k];
  avg /= double(nb);
  double var = 0.0;
  for (uint k = 1; k <= nb; k++) {
    double x = sampvals[k] - avg;
    var += x * x;
  }
  var /= double(nb - 1); // bootstrap

  const double conf_level = 0.68; // one standard deviation
  const double lowconf = 0.5 * (1.0 - conf_level);
  const double uppconf = 0.5 * (1.0 + conf_level);
  const double eps = 1e-10;
  std::sort(sampvals.c_vector_ref().begin() + 1, sampvals.c_vector_ref().end());

  // get the low and high confidence bounds and the median

  uint i = (unsigned int)floor(lowconf * nb + eps);
  double low = 0.5 * (sampvals[i] + sampvals[i + 1]);
  i = (unsigned int)floor(uppconf * nb + eps);
  double upp = 0.5 * (sampvals[i] + sampvals[i + 1]);
  i = (unsigned int)floor(0.5 * nb + eps);
  double med = 0.5 * (sampvals[i] + sampvals[i + 1]);

  result.bootassign(sampvals[0], avg, sqrt(var), low, med, upp);
}

// *************************************************************
/*

             // read all samplings from file and put into memory

void KBObsHandler::readSamplingValuesFromFile(const string& filename,
                                              XMLHandler& xmlout)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler->connectSamplingsFile(filename);
    xmlout.put_child("Status","Success");}
 catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}


             // read all samplings from file and put into memory (this version
             // only reads those records matching the MCObsInfo objects in
"obskeys")

void KBObsHandler::readSamplingValuesFromFile(const set<MCObsInfo>& obskeys,
                                              const string& filename,
                                              XMLHandler& xmlout)
{
 xmlout.set_root("ReadSamplingsFromFile");
 string fname=tidyString(filename);
 if (fname.empty()){
    xmlout.put_child("Error","Empty file name");
    return;}
 try{
    m_in_handler->connectSamplingsFile(filename,obskeys);
    xmlout.put_child("Status","Success");
    for (set<MCObsInfo>::const_iterator
it=obskeys.begin();it!=obskeys.end();it++){ XMLHandler xmlobs;
it->output(xmlobs); xmlout.put_child(xmlobs);}} catch(std::exception& xp){
    xmlout.put_child("Error",xp.what());}
}
*/
// Write from memory into file (only if all samplings available).
// If "filename" does not exist, it will be created.  Existing files are never
// erased.  New records are added to the files.  If the key of a record to be
// put already exists in a file, the put will only occur if "overwrite" is
// specified AND the size of the data to be put does not exceed the size of the
// data already in the file for that key. All observables in "obskeys" must be
// from the same ensemble.

void KBObsHandler::writeSamplingValuesToFile(const set<KBObsInfo>& obskeys,
                                             const string& filename,
                                             XMLHandler& xmlout,
                                             bool overwrite) {
  if (obskeys.empty())
    return;
  MCEnsembleInfo mcens(obskeys.begin()->getMCEnsembleInfo());
  for (set<KBObsInfo>::const_iterator it = obskeys.begin(); it != obskeys.end();
       ++it)
    if (it->getMCEnsembleInfo() != mcens)
      throw(std::runtime_error("All observables must be from same ensemble in "
                               "writeSamplingValuesToFile"));
  xmlout.set_root("WriteSamplingsToFile");
  string fname = tidyString(filename);
  if (fname.empty()) {
    xmlout.put_child("Error", "Empty file name");
    return;
  }
  xmlout.put_child("FileName", fname);
  try {
    SamplingsPutHandler SP(m_in_handler->getBinsInfo(mcens),
                           m_in_handler->getSamplingInfo(), filename,
                           overwrite);
    for (set<KBObsInfo>::const_iterator it = obskeys.begin();
         it != obskeys.end(); it++) {
      const MCObsInfo& mcobs = it->getMCObsInfo();
      XMLHandler xmlo;
      mcobs.output(xmlo);
      XMLHandler xmle("Write");
      xmle.put_child(xmlo);
      if (!queryFullAndSamplings(*it)) {
        xmle.put_child("Error", "Not all samplings are in memory");
      } else {
        try {
          const RVector& buffer = getFullAndSamplingValues(*it);
          SP.putData(mcobs, buffer);
          xmle.put_child("Success");
        } catch (std::exception& xp) {
          xmle.put_child("Error", string(xp.what()));
          xmle.put_child(xmlo);
        }
      }
      xmlout.put_child(xmle);
    }
  } catch (const std::exception& errmsg) {
    xmlout.put_child("Error", string(errmsg.what()));
  }
}

// ************************************************************************
