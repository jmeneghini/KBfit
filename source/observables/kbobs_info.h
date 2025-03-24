#ifndef KBOBS_INFO_H
#define KBOBS_INFO_H

#include "ensemble_info.h"
#include "mcobs_info.h"
#include "xml_handler.h"
#include <map>
#include <vector>

// ********************************************************************
// *                                                                  *
// *   Objects of class "KBObsInfo" store identifying information     *
// *   about one particular observable involved in KB fitting.  Each  *
// *   consists of an MCObsInfo and an MCEnsembleInfo (stored as an   *
// *   integer index).  For observables that cannot be associated     *
// *   with a single ensemble, such as K-matrix fix parameters, the   *
// *   ensemble is set to "indep" and the index is set to zero.       *
// *   Construction by XML content requires XML in the                *
// *   following format:                                              *
// *                                                                  *
// *   <KBObservable>                                                 *
// *      <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo> *
// *      <MCObservable>                                              *
// *       <ObsName>T1up_Energy</ObsName> (32 char or less, no blanks)*
// *       <Index>3</Index>        (opt nonneg integer: default 0)    *
// *       <Simple/>      (optional: if simple observable)            *
// *       <Arg>RealPart</Arg> or <Arg>Re</Arg>                       *
// *           or <Arg>ImaginaryPart</Arg> or <Arg>Im</Arg>           *
// *     </MCObservable>                                              *
// *   <KBObservable>                                                 *
// *                                                                  *
// *   If the <MCEnsembleInfo> tag is omitted, then                   *
// *          <MCEnsembleInfo>indep</MCEnsembleInfo>                  *
// *   is assumed.                                                    *
// *                                                                  *
// *   Of course, the short form of the MCObsInfo is allowed:         *
// *                                                                  *
// *   <KBObservable>                                                 *
// *      <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo> *
// *      <MCObs>T1up_Energy 3 s re</MCObs>                           *
// *   <KBObservable>                                                 *
// *                                                                  *
// ********************************************************************

// ********************************************************************
// *                                                                  *
// *   Implementation notes:                                          *
// *                                                                  *
// *   The encoding of a KBObsInfo object is nearly identical to      *
// *   that of an MCObsInfo, except the ensemble corresponding to     *
// *   each observable must be stored.  The rightmost 3 bits of       *
// *   icode[0] are kept the same as for an MCObsInfo, but the        *
// *   remaining 29 bits are changed: the first 14 bits contain an    *
// *   integer code identifying the ensemble, and the remaining       *
// *   15 bits are the MCObsInfo observable index.  A static map and  *
// *   vector are used to associate an integer index with each        *
// *   ensemble.                                                      *
// *                                                                  *
// ********************************************************************

class KBObsInfo {

  std::vector<unsigned int> icode;

  KBObsInfo(); // no default

public:
  KBObsInfo(XMLHandler& xml_in);

  KBObsInfo(const MCEnsembleInfo& mcens, const MCObsInfo& mcobs);

  KBObsInfo(const KBObsInfo& KB) : icode(KB.icode) {}

  KBObsInfo& operator=(const KBObsInfo& KB) {
    icode = KB.icode;
    return *this;
  }

  ~KBObsInfo() {}

  void setToRealPart();

  void setToImaginaryPart();

  void resetObsIndex(uint ind);

  void resetMCEnsembleInfo(const MCEnsembleInfo& mcens);

  // output functions

  const MCEnsembleInfo& getMCEnsembleInfo() const;

  MCObsInfo getMCObsInfo() const;

  std::string output(bool longform = false, int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout, bool longform = false) const; // XML output

  bool operator==(const KBObsInfo& rhs) const;

  bool operator!=(const KBObsInfo& rhs) const;

  bool operator<(const KBObsInfo& rhs) const;

private:
  static std::map<MCEnsembleInfo, uint> m_ens_map;

  static std::vector<MCEnsembleInfo> m_ensembles;
};

// ***************************************************************
#endif
