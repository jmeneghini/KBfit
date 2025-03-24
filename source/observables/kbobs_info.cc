#include "kbobs_info.h"
#include "multi_compare.h"

using namespace std;

// ***************************************************************

KBObsInfo::KBObsInfo(XMLHandler& xml_in) {
  try {
    XMLHandler xmlr(xml_in);
    xmlr.set_exceptions_on();
    XMLHandler xmlb(xmlr, "KBObservable");
    string ensid = "indep";
    xmlreadif(xmlb, "MCEnsembleInfo", ensid, "KBObsInfo");
    MCEnsembleInfo mcens(ensid);
    MCObsInfo mcobs(xmlb);
    icode = mcobs.icode;
    resetMCEnsembleInfo(mcens);
  } catch (const std::exception& msg) {
    throw(std::invalid_argument(
        string("Invalid XML for KBObsInfo constructor: ") +
        string(msg.what())));
  }
}

KBObsInfo::KBObsInfo(const MCEnsembleInfo& mcens, const MCObsInfo& mcobs)
    : icode(mcobs.icode) {
  resetMCEnsembleInfo(mcens);
}

void KBObsInfo::setToRealPart() {
  icode[0] &= ~1u; // clear the bit
}

void KBObsInfo::setToImaginaryPart() {
  icode[0] |= 1u; // set the bit
}

void KBObsInfo::resetObsIndex(uint index) {
  if (index > 32767)
    throw(std::invalid_argument("Index too large in KBObsInfo"));
  icode[0] &= ~262136u; // clear the 15 bits
  icode[0] |= (index << 3);
}

void KBObsInfo::resetMCEnsembleInfo(const MCEnsembleInfo& mcens) {
  uint ecode = 0;
  map<MCEnsembleInfo, uint>::iterator it = m_ens_map.find(mcens);
  if (it != m_ens_map.end()) {
    ecode = it->second;
  } else {
    ecode = m_ensembles.size();
    m_ensembles.push_back(mcens);
    m_ens_map.insert(make_pair(mcens, ecode));
  }
  icode[0] &= ~4294705152u;
  icode[0] |= (ecode << 18);
}

const MCEnsembleInfo& KBObsInfo::getMCEnsembleInfo() const {
  return m_ensembles[icode[0] >> 18];
}

MCObsInfo KBObsInfo::getMCObsInfo() const {
  MCObsInfo res;
  res.icode = icode;
  res.icode[0] &= ~4294705152u;
  return res;
}

string KBObsInfo::output(bool longform, int indent) const {
  XMLHandler xmlout;
  output(xmlout, longform);
  return xmlout.output(indent);
}

string KBObsInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KBObsInfo::output(XMLHandler& xmlout, bool longform) const {
  xmlout.set_root("KBObservable");
  XMLHandler xmlt;
  getMCEnsembleInfo().output(xmlt);
  xmlout.put_child(xmlt);
  getMCObsInfo().output(xmlt, longform);
  xmlout.put_child(xmlt);
}

bool KBObsInfo::operator==(const KBObsInfo& rhs) const {
  return multiEqual(icode, rhs.icode);
}

bool KBObsInfo::operator!=(const KBObsInfo& rhs) const {
  return multiNotEqual(icode, rhs.icode);
}

bool KBObsInfo::operator<(const KBObsInfo& rhs) const {
  return multiLessThan(icode, rhs.icode);
}

std::map<MCEnsembleInfo, uint> KBObsInfo::m_ens_map;
//  = {make_pair(MCEnsembleInfo("indep"),0)};

std::vector<MCEnsembleInfo> KBObsInfo::m_ensembles;
//  = {MCEnsembleInfo("indep")};

// ******************************************************************************
