#ifndef K_MATRIX_INFO_H
#define K_MATRIX_INFO_H

#include "xml_handler.h"
#include "param_registry.h"

//  This class is a compound index containing (L, 2S, a)
//     m_store encodes the L (4 bits), total S times two (4 bits),
//     channel index (4 bits)
//  Construction by XML:
//    <KIndex><L>3</L><Sx2>1<Sx2><Chan>0</Chan></KIndex>
//  or <KIndex>L(3) 2S(1) chan(0)</KIndex>  (order matters in this form)

class KIndex {
  uint m_store;

  KIndex() = delete; // no default value

public:
  KIndex(XMLHandler& xmlin);

  KIndex(uint L, uint Stimestwo, uint channel_index);

  KIndex(const KIndex& in) : m_store(in.m_store) {}

  KIndex& operator=(const KIndex& in) {
    m_store = in.m_store;
    return *this;
  }

  ~KIndex() {}

  uint getL() const;

  uint getStimestwo() const;

  uint getChannelIndex() const;

  bool operator==(const KIndex& rhs) const;

  bool operator!=(const KIndex& rhs) const;

  bool operator<(const KIndex& rhs) const;

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

private:
  KIndex(uint store) : m_store(store) {}

  void check_encode(uint L, uint Stimestwo, uint channel_index);

  friend class KElementInfo;

  friend class KFitParamInfo;
};

inline uint KIndex::getL() const { return (m_store >> 8) & 0xFu; }

inline uint KIndex::getStimestwo() const { return (m_store >> 4) & 0xFu; }

inline uint KIndex::getChannelIndex() const { return m_store & 0xFu; }

inline bool KIndex::operator==(const KIndex& rhs) const {
  return (m_store == rhs.m_store);
}

inline bool KIndex::operator!=(const KIndex& rhs) const {
  return (m_store != rhs.m_store);
}

inline bool KIndex::operator<(const KIndex& rhs) const {
  return (m_store < rhs.m_store);
}

//  This class identifies an element of the K-matrix using
//          J (L'S'a', LSa)
//     m_store encodes
//       J (four bits)  (L'S'a') (12 bits),  (Lsa) (12 bits)
//  Since K is real and symmetric, the (L'S'a') and (LSa)
//  order is sorted according to how stored in "KIndex".
//  The constructor also checks that J,L,S are allowed values.
//  XML input for constructor must have form
//    <KElementInfo>
//       <JTimesTwo>6</JTimesTwo>
//       <KIndex>...</KIndex>
//       <KIndex>...</KIndex>    ( must have 1 or 2 <KIndex> tags)
//    </KElementInfo>

class KElementInfo {
  uint m_store;

  KElementInfo() = delete; // no default value

public:
  KElementInfo(XMLHandler& xmlin);

  KElementInfo(uint Jtimestwo, uint rowL, uint rowStimestwo,
               uint rowchannelindex, uint colL, uint colStimestwo,
               uint colchannelindex);

  KElementInfo(uint Jtimestwo, const KIndex& row, const KIndex& col);

  KElementInfo(const KElementInfo& in) : m_store(in.m_store) {}

  KElementInfo& operator=(const KElementInfo& in) {
    m_store = in.m_store;
    return *this;
  }

  ~KElementInfo() {}

  uint getJtimestwo() const;

  uint getRowL() const;

  uint getRowStimestwo() const;

  uint getRowChannelIndex() const;

  KIndex getRow() const;

  uint getColumnL() const;

  uint getColumnStimestwo() const;

  uint getColumnChannelIndex() const;

  KIndex getColumn() const;

  bool operator==(const KElementInfo& rhs) const;

  bool operator!=(const KElementInfo& rhs) const;

  bool operator<(const KElementInfo& rhs) const;

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

private:
  void check_encode(uint Jtimestwo, uint rowL, uint rowStimestwo,
                    uint rowchannelindex, uint colL, uint colStimestwo,
                    uint colchannelindex);

  void encode(uint Jtimestwo, const KIndex& kl, const KIndex& kr);

  bool invalidJLS(uint Jtimestwo, uint L, uint Stimestwo);

  KElementInfo(uint store) : m_store(store) {}

  friend class KFitParamInfo;
};

inline uint KElementInfo::getJtimestwo() const {
  return (m_store >> 24) & 0xFu;
}

inline uint KElementInfo::getRowL() const { return (m_store >> 20) & 0xFu; }

inline uint KElementInfo::getRowStimestwo() const {
  return (m_store >> 16) & 0xFu;
}

inline uint KElementInfo::getRowChannelIndex() const {
  return (m_store >> 12) & 0xFu;
}

inline KIndex KElementInfo::getRow() const {
  return KIndex((m_store >> 12) & 0xFFFu);
}

inline uint KElementInfo::getColumnL() const { return (m_store >> 8) & 0xFu; }

inline uint KElementInfo::getColumnStimestwo() const {
  return (m_store >> 4) & 0xFu;
}

inline uint KElementInfo::getColumnChannelIndex() const {
  return m_store & 0xFu;
}

inline KIndex KElementInfo::getColumn() const {
  return KIndex(m_store & 0xFFFu);
}

inline bool KElementInfo::operator==(const KElementInfo& rhs) const {
  return (m_store == rhs.m_store);
}

inline bool KElementInfo::operator!=(const KElementInfo& rhs) const {
  return (m_store != rhs.m_store);
}

inline bool KElementInfo::operator<(const KElementInfo& rhs) const {
  return (m_store < rhs.m_store);
}

//  This stores information about a given K-matrix fit parameter.
//  This class is mainly informational for the end user.
//
//    m_store1: left-most 4 bits encode the type of parameter
//         0 = term in a polynomial
//         1 = pole energy
//         2 = coupling in pole residue
//
//    if polynomial:
//         - rightmost 28 bits of m_store1 contains "n"
//               where parameter is coefficient of Ecm^n
//         - m_store2 contains the m_store of the K-matrix element
//               the polynomial is associated with
//
//    if pole energy:
//         - rightmost 28 bits of m_store1 contains the integer index
//            identifying the pole (14 bits) and total J times two (14 bits)
//         - m_store2 = 0
//
//    if coupling in pole residue
//         - rightmost 28 bits of m_store1 contain integer index
//              identify the pole (14 bits) and total J times two (14 bits)
//         - m_store2 contains m_store of the KIndex associated
//              with this parameter
//
//  Construction by XML input of the form:
//
//    <KFitParamInfo>
//      <PolynomialTerm>
//         <Power>3</Power>
//         <KElementInfo>...</KElementInfo>
//      </PolynomialTerm>
//    </KFitParamInfo>
//
//    <KFitParamInfo>
//      <PoleEnergy>
//         <Index>3</Index>
//         <JTimesTwo>2</JTimesTwo>
//      </PoleEnergy>
//    </KFitParamInfo>
//
//    <KFitParamInfo>
//      <PoleCoupling>
//         <Index>3</Index>
//         <JTimesTwo>2</JTimesTwo>
//         <KIndex>...</KIndex>
//      </PoleCoupling>
//    </KFitParamInfo>
//
//  Two routines that are needed for converting between a KFitParamInfo
//  object and an MCObsInfo (as needed for a record key in a samplings file)
//  are
//
//      void setFromMCObsName(const std::string& obsname);
//      std::string getMCObsName() const;
//
//  To make an MCObsInfo from a KFitParamInfo, a string of 32 or less
//  characters is needed for the MCObsInfo name.  These are chosen here
//  as
//
//     KPoly(n,2J)(L',2S',a')(L,2S,a)   --> polynomial term  n = power
//
//     KPoleE(n,2J)   --> pole energy   n = index of pole
//
//     KPoleC(n,2J)(L,2S,a)     --> pole coupling  n = index of pole
//

class KFitParamInfo {

  uint m_store1;
  uint m_store2;

public:
  KFitParamInfo(); // default 0,0

  KFitParamInfo(XMLHandler& xmlin);

  KFitParamInfo(const KElementInfo& keleminfo, uint power);

  KFitParamInfo(uint pole_index, uint Jtimestwo);

  KFitParamInfo(const KIndex& kindex, uint pole_index, uint Jtimestwo);

  // using hash instead of uint to have a constructor for expressions that
  // doesn't conflict with the above polynomial constructor
  KFitParamInfo(const KElementInfo& keleminfo, std::string param_name,
                ParameterNameRegistry& param_registry);

  KFitParamInfo(const KFitParamInfo& in);

  KFitParamInfo& operator=(const KFitParamInfo& in);

  ~KFitParamInfo() {}

  void setPolynomialPower(const KElementInfo& keleminfo, uint power);

  void setPoleEnergy(uint pole_index, uint Jtimestwo);

  void setPoleCoupling(const KIndex& kindex, uint pole_index, uint Jtimestwo);

  void setFromMCObsName(const std::string& obsname);

  bool operator==(const KFitParamInfo& rhs) const;

  bool operator!=(const KFitParamInfo& rhs) const;

  bool operator<(const KFitParamInfo& rhs) const;

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

  std::string getMCObsName() const;

private:
  void set_poly_power(const KElementInfo& keleminfo, uint power);

  void set_pole_energy(uint pole_index, uint Jtimestwo);

  void set_pole_coupling(const KIndex& kindex, uint pole_index, uint Jtimestwo);

  void set_string_expr_param(const KElementInfo& keleminfo, uint param_hash);
};

inline KFitParamInfo::KFitParamInfo()
    : m_store1(0), m_store2(0) {} // default value

inline KFitParamInfo::KFitParamInfo(const KElementInfo& keleminfo, uint power) {
  set_poly_power(keleminfo, power);
}

inline KFitParamInfo::KFitParamInfo(uint pole_index, uint Jtimestwo) {
  set_pole_energy(pole_index, Jtimestwo);
}

inline KFitParamInfo::KFitParamInfo(const KIndex& kindex, uint pole_index,
                                    uint Jtimestwo) {
  set_pole_coupling(kindex, pole_index, Jtimestwo);
}

inline KFitParamInfo::KFitParamInfo(const KElementInfo& keleminfo, std::string param_name,
                                    ParameterNameRegistry& param_registry) {
  uint param_hash = param_registry.registerParameter(param_name);
  set_string_expr_param(keleminfo, param_hash);
}

inline KFitParamInfo::KFitParamInfo(const KFitParamInfo& in)
    : m_store1(in.m_store1), m_store2(in.m_store2) {}

inline KFitParamInfo& KFitParamInfo::operator=(const KFitParamInfo& in) {
  m_store1 = in.m_store1;
  m_store2 = in.m_store2;
  return *this;
}

inline void KFitParamInfo::setPolynomialPower(const KElementInfo& keleminfo,
                                              uint power) {
  set_poly_power(keleminfo, power);
}

inline void KFitParamInfo::setPoleEnergy(uint pole_index, uint Jtimestwo) {
  set_pole_energy(pole_index, Jtimestwo);
}

inline void KFitParamInfo::setPoleCoupling(const KIndex& kindex,
                                           uint pole_index, uint Jtimestwo) {
  set_pole_coupling(kindex, pole_index, Jtimestwo);
}

inline bool KFitParamInfo::operator==(const KFitParamInfo& rhs) const {
  return (m_store1 == rhs.m_store1) && (m_store2 == rhs.m_store2);
}

inline bool KFitParamInfo::operator!=(const KFitParamInfo& rhs) const {
  return (m_store1 != rhs.m_store1) || (m_store2 != rhs.m_store2);
}

inline bool KFitParamInfo::operator<(const KFitParamInfo& rhs) const {
  return (m_store1 < rhs.m_store1) ||
         ((m_store1 == rhs.m_store1) && (m_store2 < rhs.m_store2));
}

// *************************************************************************

//  Contains information related to a particular decay channel:
//  names of each particle (such as "pion", "eta", "pion+"),
//  spins of each particle (times two to make a nonnegative integer),
//  and booleans telling if the particles are identical and if they
//  have the same intrinsic parities.

//  Construction by XML:
//    <DecayChannelInfo>
//      <Particle1Name>pion</Particle1Name>
//      <Spin1TimesTwo>0</Spin1TimesTwo>
//      <Identical/> (if identical, do not include tags below)
//      <Particle2Name>eta</Particle2Name>
//      <Spin2TimesTwo>2</Spin2TimesTwo>
//      <IntrinsicParities>same</IntrinsicParities> (or "opposite")
//    </DecayChannelInfo>

//  Implementation: m_store contains 15 bits spin1timestwo,
//  15 bits spin2timestwo, 1 bit identical, 1 bit same intrinsic parities.

class DecayChannelInfo {
  std::string m_name1, m_name2;
  uint m_store;

public:
  DecayChannelInfo(const std::string& particle1_name,
                   const std::string& particle2_name, uint spin1_times_two,
                   uint spin2_times_two, bool identical,
                   bool same_intrinsic_parities);

  DecayChannelInfo(XMLHandler& xmlin);

  DecayChannelInfo(const DecayChannelInfo& dcinfo);

  DecayChannelInfo& operator=(const DecayChannelInfo& dcinfo);

  ~DecayChannelInfo() {}

  std::string getParticle1Name() const { return m_name1; }

  std::string getParticle2Name() const { return m_name2; }

  uint getSpin1timestwo() const;

  uint getSpin2timestwo() const;

  bool areIdentical() const;

  bool sameIntrinsicParities() const;

  std::string output(int indent = 0) const; // XML output

  std::string str() const; // XML output

  void output(XMLHandler& xmlout) const; // XML output

  bool operator==(const DecayChannelInfo& rhs) const;

  bool operator!=(const DecayChannelInfo& rhs) const;

  // bool operator<(const DecayChannelInfo& rhs) const;

private:
  void check_encode(uint spin1_times_two, uint spin2_times_two, bool identical,
                    bool same_intrinsic_parities);
};

std::ostream& operator<<(std::ostream& os, const DecayChannelInfo& dcinfo);

inline uint DecayChannelInfo::getSpin1timestwo() const { return m_store >> 17; }

inline uint DecayChannelInfo::getSpin2timestwo() const {
  return (m_store >> 2) & 0x7FFFu;
}

inline bool DecayChannelInfo::areIdentical() const {
  return (m_store >> 1) & 0x1u;
}

inline bool DecayChannelInfo::sameIntrinsicParities() const {
  return m_store & 0x1u;
}

inline bool DecayChannelInfo::operator==(const DecayChannelInfo& rhs) const {
  return (m_store == rhs.m_store) && (m_name1 == rhs.m_name1) &&
         (m_name2 == rhs.m_name2) && (m_name1 == rhs.m_name1) &&
         (m_name2 == rhs.m_name2);
}

inline bool DecayChannelInfo::operator!=(const DecayChannelInfo& rhs) const {
  return (m_store != rhs.m_store) || (m_name1 != rhs.m_name1) ||
         (m_name2 != rhs.m_name2) || (m_name1 != rhs.m_name1) ||
         (m_name2 != rhs.m_name2);
}

// *************************************************************************
#endif
