#include "K_matrix_calc.h"

using namespace std;

// ****************************************************************

KtildeMatrixCalculator::~KtildeMatrixCalculator() {
  for (std::map<KElementInfo, FitForm*>::iterator it = m_fit.begin();
       it != m_fit.end(); it++)
    delete it->second;
}

KtildeMatrixCalculator::KtildeMatrixCalculator(XMLHandler& xmlin,
                                               bool require_initvals) {
  XMLHandler xmlk(xmlin, "KtildeMatrix");
  XMLHandler xmld(xmlk, "DecayChannels");
  list<XMLHandler> dcxml = xmld.find_among_children("DecayChannelInfo");
  for (list<XMLHandler>::iterator it = dcxml.begin(); it != dcxml.end(); it++) {
    DecayChannelInfo dctemp(*it);
    m_decay_infos.push_back(dctemp);
  }
  if (m_decay_infos.empty())
    throw(std::invalid_argument("Must have at least ONE decay channel"));
  list<XMLHandler> kelems = xmlk.find_among_children("Element");
  initialize(kelems);
  if (require_initvals)
    initialize_starting_values(xmlk);
}

uint KtildeMatrixCalculator::getNumberOfParameters() const {
  return m_paraminfo.size();
}

void KtildeMatrixCalculator::setParameterValues(
    std::vector<double> kappa_params) {
  if (kappa_params.size() != getNumberOfParameters())
    throw(std::invalid_argument("Could not set KtildeParameters"));
  m_kappa_params = kappa_params;
}

int KtildeMatrixCalculator::getParameterIndex(
    const KFitParamInfo& kinfo) const {
  std::map<KFitParamInfo, uint>::const_iterator it;
  it = m_paramindices.find(kinfo);
  if (it == m_paramindices.end())
    return -1;
  return it->second;
}

double
KtildeMatrixCalculator::getParameterValue(const KFitParamInfo& kinfo) const {
  return m_kappa_params.at(getParameterIndex(kinfo));
}

set<KElementInfo> KtildeMatrixCalculator::getElementInfos() const {
  set<KElementInfo> res;
  std::map<KElementInfo, FitForm*>::const_iterator it;
  for (it = m_fit.begin(); it != m_fit.end(); it++)
    res.insert(it->first);
  return res;
}

string KtildeMatrixCalculator::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string KtildeMatrixCalculator::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KtildeMatrixCalculator::output(XMLHandler& xmlout) const {
  xmlout.set_root("KtildeMatrixCalculator");
  XMLHandler xmlk("KtildeMatrix");
  XMLHandler xmld("DecayChannels");
  for (uint chan = 0; chan < m_decay_infos.size(); chan++) {
    XMLHandler xmlci;
    m_decay_infos[chan].output(xmlci);
    xmld.put_child(xmlci);
  }
  xmlk.put_child(xmld);
  for (map<KElementInfo, FitForm*>::const_iterator it = m_fit.begin();
       it != m_fit.end(); it++) {
    XMLHandler xmle;
    (it->first).output(xmle);
    XMLHandler xmlf;
    (it->second)->output(xmlf);
    XMLHandler xmlc("Element");
    xmlc.put_child(xmle);
    xmlc.put_child(xmlf);
    xmlk.put_child(xmlc);
  }
  xmlout.put_child(xmlk);
  XMLHandler xmlps("FitParameterInfos");
  for (uint k = 0; k < m_paraminfo.size(); k++) {
    XMLHandler xmlp("FitParameterInfo");
    XMLHandler xmlpp;
    m_paraminfo[k].output(xmlpp);
    xmlp.put_child(xmlpp);
    xmlp.put_child("Index", make_string(k));
    xmlps.put_child(xmlp);
  }
  xmlout.put_child(xmlps);
}

void KtildeMatrixCalculator::initialize(list<XMLHandler>& kelems) {
  try {
    uint nchan = m_decay_infos.size();
    std::map<KFitParamInfo, double> initvals;
    for (list<XMLHandler>::iterator it = kelems.begin(); it != kelems.end();
         it++) {
      KElementInfo keinfo(*it);
      if ((keinfo.getRowChannelIndex() >= nchan) ||
          (keinfo.getColumnChannelIndex() >= nchan))
        throw(
            std::invalid_argument("KElementInfo has invalid channel indices"));
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      initialize_a_fitform(keinfo, *it);
    }
    m_paraminfo.resize(m_paramindices.size());
    for (std::map<KFitParamInfo, uint>::iterator it = m_paramindices.begin();
         it != m_paramindices.end(); it++)
      m_paraminfo[it->second] = it->first;
    m_kappa_params.resize(m_paraminfo.size());
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not construct Kmatrix: ") +
                                xp.what()));
  }
}

//  very important class that sets up m_fit and m_paramindices

void KtildeMatrixCalculator::initialize_a_fitform(const KElementInfo& kinfo,
                                                  XMLHandler& xmlin) {
  FitForm* fptr = 0;
  try {
    XMLHandler xmle(xmlin, "FitForm");
    uint count1 = xmle.count_among_children("Polynomial");
    uint count2 = xmle.count_among_children("SumOfPoles");
    uint count3 = xmle.count_among_children("SumOfPolesPlusPolynomial");
    uint count4 = xmle.count_among_children("Expression");
    if ((count1 + count2 + count3 + count4) != 1)
      throw(std::invalid_argument(string("FitForm has invalid XML content: ") +
                                  xmle.str()));
    if (count1 == 1) {
      XMLHandler xmlf(xmle, "Polynomial");
      fptr = new Polynomial(xmlf);
    } else if (count2 == 1) {
      XMLHandler xmlf(xmle, "SumOfPoles");
      fptr = new SumOfPoles(xmlf);
    } else if (count3 == 1) {
      XMLHandler xmlf(xmle, "SumOfPolesPlusPolynomial");
      fptr = new SumOfPolesPlusPolynomial(xmlf);
    } else {
      XMLHandler xmlf(xmle, "Expression");
      fptr = new Expression(xmlf);
    }
    fptr->Kinitialize(kinfo, m_paramindices);
    m_fit.insert(make_pair(kinfo, fptr));
  }

  catch (const std::exception& xp) {
    delete fptr;
    throw(std::invalid_argument(
        string("In KtildeMatrix, could not initialize a fit form: ") +
        xp.what()));
  }
}

void KtildeMatrixCalculator::initialize_starting_values(XMLHandler& xmlin) {
  try {
    double val = 0.0;
    XMLHandler xmls(xmlin, "StartingValues");
    list<XMLHandler> xmlp = xmls.find("KFitParamInfo");
    map<KFitParamInfo, double> initvals;
    for (list<XMLHandler>::iterator it = xmlp.begin(); it != xmlp.end(); ++it) {
      KFitParamInfo kpinfo(*it);
      xmlread(*it, "StartingValue", val, "KtildeMatrixCalculator");
      initvals[kpinfo] = val;
    }
    for (uint k = 0; k < m_paraminfo.size(); ++k) {
      std::map<KFitParamInfo, double>::const_iterator kt =
          initvals.find(m_paraminfo[k]);
      if (kt == initvals.end())
        throw(
            std::invalid_argument(string("Could not find initial value for ") +
                                  m_paraminfo[k].str()));
      m_kappa_params[k] = kt->second;
    }
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(
        string("Could not initialize starting values: ") + xp.what()));
  }
}

KtildeMatrixCalculator::KtildeMatrixCalculator(
    const std::list<std::pair<KElementInfo, Polynomial>>& pelems,
    const std::list<std::pair<KElementInfo, SumOfPoles>>& selems,
    const std::list<std::pair<KElementInfo, SumOfPolesPlusPolynomial>>& spelems,
    const std::vector<DecayChannelInfo>& chans) {
  m_decay_infos = chans;
  FitForm* fptr = 0;
  try {
    for (std::list<std::pair<KElementInfo, Polynomial>>::const_iterator it =
             pelems.begin();
         it != pelems.end(); it++) {
      const KElementInfo& keinfo(it->first);
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      fptr = new Polynomial(it->second);
      fptr->Kinitialize(keinfo, m_paramindices);
      m_fit.insert(make_pair(keinfo, fptr));
    }
    for (std::list<std::pair<KElementInfo, SumOfPoles>>::const_iterator it =
             selems.begin();
         it != selems.end(); it++) {
      const KElementInfo& keinfo(it->first);
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      fptr = new SumOfPoles(it->second);
      fptr->Kinitialize(keinfo, m_paramindices);
      m_fit.insert(make_pair(keinfo, fptr));
    }
    for (std::list<std::pair<KElementInfo, SumOfPolesPlusPolynomial>>::
             const_iterator it = spelems.begin();
         it != spelems.end(); it++) {
      const KElementInfo& keinfo(it->first);
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      fptr = new SumOfPolesPlusPolynomial(it->second);
      fptr->Kinitialize(keinfo, m_paramindices);
      m_fit.insert(make_pair(keinfo, fptr));
    }
    m_paraminfo.resize(m_paramindices.size());
    for (std::map<KFitParamInfo, uint>::iterator it = m_paramindices.begin();
         it != m_paramindices.end(); it++)
      m_paraminfo[it->second] = it->first;
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not construct Kmatrix: ") +
                                xp.what()));
  }
}

double KtildeMatrixCalculator::calculate(uint Jtimestwo, uint Lp,
                                         uint Sptimestwo, uint chanp, uint L,
                                         uint Stimestwo, uint chan,
                                         double Ecm_over_mref) const {
  KElementInfo kelem(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan);
  std::map<KElementInfo, FitForm*>::const_iterator it = m_fit.find(kelem);
  if (it == m_fit.end())
    return 0.0;
  return (it->second)->evaluate(m_kappa_params, Ecm_over_mref);
}

bool KtildeMatrixCalculator::isZero(uint Jtimestwo, uint Lp, uint Sptimestwo,
                                    uint chanp, uint L, uint Stimestwo,
                                    uint chan) const {
  KElementInfo kelem(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan);
  std::map<KElementInfo, FitForm*>::const_iterator it = m_fit.find(kelem);
  return (it == m_fit.end());
}

// ***************************************************************************************

KtildeInverseCalculator::~KtildeInverseCalculator() {
  for (std::map<KElementInfo, FitForm*>::iterator it = m_fit.begin();
       it != m_fit.end(); it++)
    delete it->second;
}

KtildeInverseCalculator::KtildeInverseCalculator(XMLHandler& xmlin,
                                                 bool require_initvals) {
  XMLHandler xmlk(xmlin, "KtildeMatrixInverse");
  XMLHandler xmld(xmlk, "DecayChannels");
  list<XMLHandler> dcxml = xmld.find_among_children("DecayChannelInfo");
  for (list<XMLHandler>::iterator it = dcxml.begin(); it != dcxml.end(); it++) {
    DecayChannelInfo dctemp(*it);
    m_decay_infos.push_back(dctemp);
  }
  if (m_decay_infos.empty())
    throw(std::invalid_argument("Must have at least ONE decay channel"));
  list<XMLHandler> kelems = xmlk.find_among_children("Element");
  initialize(kelems);
  if (require_initvals)
    initialize_starting_values(xmlk);
}

uint KtildeInverseCalculator::getNumberOfParameters() const {
  return m_paraminfo.size();
}

void KtildeInverseCalculator::setParameterValues(
    std::vector<double> kappa_params) {
  if (kappa_params.size() != getNumberOfParameters())
    throw(std::invalid_argument("Could not set KtildeParameters"));
  m_kappa_params = kappa_params;
}

int KtildeInverseCalculator::getParameterIndex(
    const KFitParamInfo& kinfo) const {
  std::map<KFitParamInfo, uint>::const_iterator it;
  it = m_paramindices.find(kinfo);
  if (it == m_paramindices.end())
    return -1;
  return it->second;
}

double
KtildeInverseCalculator::getParameterValue(const KFitParamInfo& kinfo) const {
  return m_kappa_params.at(getParameterIndex(kinfo));
}

set<KElementInfo> KtildeInverseCalculator::getElementInfos() const {
  set<KElementInfo> res;
  std::map<KElementInfo, FitForm*>::const_iterator it;
  for (it = m_fit.begin(); it != m_fit.end(); it++)
    res.insert(it->first);
  return res;
}

string KtildeInverseCalculator::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string KtildeInverseCalculator::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KtildeInverseCalculator::output(XMLHandler& xmlout) const {
  xmlout.set_root("KtildeInverseCalculator");
  XMLHandler xmlk("KtildeMatrixInverse");
  XMLHandler xmld("DecayChannels");
  for (uint chan = 0; chan < m_decay_infos.size(); chan++) {
    XMLHandler xmlci;
    m_decay_infos[chan].output(xmlci);
    xmld.put_child(xmlci);
  }
  xmlk.put_child(xmld);
  for (map<KElementInfo, FitForm*>::const_iterator it = m_fit.begin();
       it != m_fit.end(); it++) {
    XMLHandler xmle;
    (it->first).output(xmle);
    XMLHandler xmlf;
    (it->second)->output(xmlf);
    XMLHandler xmlc("Element");
    xmlc.put_child(xmle);
    xmlc.put_child(xmlf);
    xmlk.put_child(xmlc);
  }
  xmlout.put_child(xmlk);
  XMLHandler xmlps("FitParameterInfos");
  for (uint k = 0; k < m_paraminfo.size(); k++) {
    XMLHandler xmlp("FitParameterInfo");
    XMLHandler xmlpp;
    m_paraminfo[k].output(xmlpp);
    xmlp.put_child(xmlpp);
    xmlp.put_child("Index", make_string(k));
    xmlps.put_child(xmlp);
  }
  xmlout.put_child(xmlps);
}

void KtildeInverseCalculator::initialize(list<XMLHandler>& kelems) {
  try {
    uint nchan = m_decay_infos.size();
    std::map<KFitParamInfo, double> initvals;
    for (list<XMLHandler>::iterator it = kelems.begin(); it != kelems.end();
         it++) {
      KElementInfo keinfo(*it);
      if ((keinfo.getRowChannelIndex() >= nchan) ||
          (keinfo.getColumnChannelIndex() >= nchan))
        throw(
            std::invalid_argument("KElementInfo has invalid channel indices"));
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      initialize_a_fitform(keinfo, *it);
    }
    m_paraminfo.resize(m_paramindices.size());
    for (std::map<KFitParamInfo, uint>::iterator it = m_paramindices.begin();
         it != m_paramindices.end(); it++)
      m_paraminfo[it->second] = it->first;
    m_kappa_params.resize(m_paraminfo.size());
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not construct Kmatrix: ") +
                                xp.what()));
  }
}

//  very important class that sets up m_fit and m_paramindices
//  written as static class so can be used by KtildeInverseCalculator too

void KtildeInverseCalculator::initialize_a_fitform(const KElementInfo& kinfo,
                                                   XMLHandler& xmlin) {
  FitForm* fptr = 0;
  try {
    XMLHandler xmle(xmlin, "FitForm");
    uint count1 = xmle.count_among_children("Polynomial");
    uint count2 = xmle.count_among_children("Expression");
    if ((count1 + count2) != 1)
      throw(std::invalid_argument(string("FitForm has invalid XML content: ") +
                                  xmle.str()));
    if (count1 == 1) {
      XMLHandler xmlf(xmle, "Polynomial");
      fptr = new Polynomial(xmlf);
    } else {
      XMLHandler xmlf(xmle, "Expression");
      fptr = new Expression(xmlf);
    }
    fptr->Kinitialize(kinfo, m_paramindices);
    m_fit.insert(make_pair(kinfo, fptr));
  }

  catch (const std::exception& xp) {
    delete fptr;
    throw(std::invalid_argument(
        string("In KtildeInverse, could not initialize a fit form: ") +
        xp.what()));
  }
}

void KtildeInverseCalculator::initialize_starting_values(XMLHandler& xmlin) {
  try {
    double val = 0.0;
    XMLHandler xmls(xmlin, "StartingValues");
    list<XMLHandler> xmlp = xmls.find("KFitParamInfo");
    map<KFitParamInfo, double> initvals;
    for (list<XMLHandler>::iterator it = xmlp.begin(); it != xmlp.end(); ++it) {
      KFitParamInfo kpinfo(*it);
      xmlread(*it, "StartingValue", val, "KtildeInverseCalculator");
      initvals[kpinfo] = val;
    }
    for (uint k = 0; k < m_paraminfo.size(); ++k) {
      std::map<KFitParamInfo, double>::const_iterator kt =
          initvals.find(m_paraminfo[k]);
      if (kt == initvals.end())
        throw(
            std::invalid_argument(string("Could not find initial value for ") +
                                  m_paraminfo[k].str()));
      m_kappa_params[k] = kt->second;
    }
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(
        string("Could not initialize starting values: ") + xp.what()));
  }
}

KtildeInverseCalculator::KtildeInverseCalculator(
    const std::list<std::pair<KElementInfo, Polynomial>>& pelems,
    const std::vector<DecayChannelInfo>& chans) {
  m_decay_infos = chans;
  FitForm* fptr = 0;
  try {
    for (std::list<std::pair<KElementInfo, Polynomial>>::const_iterator it =
             pelems.begin();
         it != pelems.end(); it++) {
      const KElementInfo& keinfo(it->first);
      if (m_fit.find(keinfo) != m_fit.end())
        throw(std::invalid_argument("Duplicate KElementInfo tag"));
      fptr = new Polynomial(it->second);
      fptr->Kinitialize(keinfo, m_paramindices);
      m_fit.insert(make_pair(keinfo, fptr));
    }
    m_paraminfo.resize(m_paramindices.size());
    for (std::map<KFitParamInfo, uint>::iterator it = m_paramindices.begin();
         it != m_paramindices.end(); it++)
      m_paraminfo[it->second] = it->first;
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not construct Kmatrix: ") +
                                xp.what()));
  }
}

double KtildeInverseCalculator::calculate(uint Jtimestwo, uint Lp,
                                          uint Sptimestwo, uint chanp, uint L,
                                          uint Stimestwo, uint chan,
                                          double Ecm_over_mref) const {
  KElementInfo kelem(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan);
  std::map<KElementInfo, FitForm*>::const_iterator it = m_fit.find(kelem);
  if (it == m_fit.end())
    return 0.0;
  return (it->second)->evaluate(m_kappa_params, Ecm_over_mref);
}

bool KtildeInverseCalculator::isZero(uint Jtimestwo, uint Lp, uint Sptimestwo,
                                     uint chanp, uint L, uint Stimestwo,
                                     uint chan) const {
  KElementInfo kelem(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan);
  std::map<KElementInfo, FitForm*>::const_iterator it = m_fit.find(kelem);
  return (it == m_fit.end());
}

// ***************************************************************************************
