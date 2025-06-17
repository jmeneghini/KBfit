#include "K_matrix_info.h"

using namespace std;

// ****************************************************************

KIndex::KIndex(uint L, uint Stimestwo, uint channel_index) {
  check_encode(L, Stimestwo, channel_index);
}

KIndex::KIndex(XMLHandler& xmlin) {
  uint L, Stimestwo, chan;
  XMLHandler xmlk(xmlin, "KIndex");
  uint count = 0;
  if (xmlreadifchild(xmlk, "L", L))
    count++;
  if (xmlreadifchild(xmlk, "Sx2", Stimestwo))
    count++;
  if (xmlreadifchild(xmlk, "Chan", chan))
    count++;
  if (count == 0) {
    string kstr;
    xmlread(xmlin, "KIndex", kstr, "KIndex");
    count =
        sscanf(kstr.c_str(), "L(%d) 2S(%d) chan(%d)", &L, &Stimestwo, &chan);
  }
  if (count != 3)
    throw(std::invalid_argument("Invalid XML for KIndex"));
  check_encode(L, Stimestwo, chan);
}

//   m_store encodes the channel index (6 bits),
//   total S times two (6 bits), total J times two (7 bits)
//   L (7 bits), occurrence index (6 bits)

void KIndex::check_encode(uint L, uint Stimestwo, uint channel_index) {
  if ((channel_index >= 16) || (Stimestwo >= 16) || (L >= 16))
    throw(std::invalid_argument("Unsupported KIndex value"));
  m_store = L;
  m_store <<= 4;
  m_store |= Stimestwo;
  m_store <<= 4;
  m_store |= channel_index;
}

string KIndex::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string KIndex::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KIndex::output(XMLHandler& xmlout) const {
  string ststr(" L(");
  ststr += make_string(getL()) + ")";
  ststr += " 2S(" + make_string(getStimestwo()) + ")";
  ststr += " chan(" + make_string(getChannelIndex()) + ")";
  xmlout.set_root("KIndex", ststr);
}

// ****************************************************************

KElementInfo::KElementInfo(uint Jtimestwo, uint rowL, uint rowStimestwo,
                           uint rowchannelindex, uint colL, uint colStimestwo,
                           uint colchannelindex) {
  check_encode(Jtimestwo, rowL, rowStimestwo, rowchannelindex, colL,
               colStimestwo, colchannelindex);
}

KElementInfo::KElementInfo(uint Jtimestwo, const KIndex& row,
                           const KIndex& col) {
  if (row < col)
    encode(Jtimestwo, row, col);
  else
    encode(Jtimestwo, col, row);
}

KElementInfo::KElementInfo(XMLHandler& xmlin) {
  XMLHandler xmlk(xmlin, "KElementInfo");
  uint Jtimestwo;
  xmlreadchild(xmlk, "JTimesTwo", Jtimestwo);
  list<XMLHandler> ki = xmlk.find_among_children("KIndex");
  if ((ki.size() < 1) || (ki.size() > 2))
    throw(std::invalid_argument("Invalid XML to KElementInfo"));
  uint rowL, rowStimestwo, rowchan, colL, colStimestwo, colchan;
  KIndex row(ki.front());
  rowL = row.getL();
  rowStimestwo = row.getStimestwo();
  rowchan = row.getChannelIndex();
  if (ki.size() == 1) {
    colL = rowL;
    colStimestwo = rowStimestwo;
    colchan = rowchan;
  } else {
    KIndex col(*(++(ki.begin())));
    colL = col.getL();
    colStimestwo = col.getStimestwo();
    colchan = col.getChannelIndex();
  }
  check_encode(Jtimestwo, rowL, rowStimestwo, rowchan, colL, colStimestwo,
               colchan);
}

void KElementInfo::check_encode(uint Jtimestwo, uint rowL, uint rowStimestwo,
                                uint rowchannelindex, uint colL,
                                uint colStimestwo, uint colchannelindex) {
  if ((invalidJLS(Jtimestwo, rowL, rowStimestwo)) or
      (invalidJLS(Jtimestwo, colL, colStimestwo)))
    throw(std::invalid_argument("Invalid JLS in KElementInfo"));
  KIndex krow(rowL, rowStimestwo, rowchannelindex);
  KIndex kcol(colL, colStimestwo, colchannelindex);
  if (krow < kcol)
    encode(Jtimestwo, krow, kcol);
  else
    encode(Jtimestwo, kcol, krow);
}

//     m_store encodes
//       J (four bits)  (L'S'a') (12 bits),  (Lsa) (12 bits)
//  Since K is real and symmetric, the (L'S'a') and (LSa)
//  order is sorted according to how stored in "KIndex"

void KElementInfo::encode(uint Jtimestwo, const KIndex& kl, const KIndex& kr) {
  if (Jtimestwo >= 16)
    throw(std::invalid_argument("Unsupported KElementInfo value"));
  m_store = Jtimestwo;
  m_store <<= 12;
  m_store |= kl.m_store;
  m_store <<= 12;
  m_store |= kr.m_store;
}

bool KElementInfo::invalidJLS(uint Jtimestwo, uint L, uint Stimestwo) {
  return (((Jtimestwo % 2) != (Stimestwo % 2)) ||
          ((2 * L) < ((Jtimestwo >= Stimestwo) ? (Jtimestwo - Stimestwo)
                                               : (Stimestwo - Jtimestwo))) ||
          ((2 * L) > (Jtimestwo + Stimestwo)));
}

string KElementInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string KElementInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KElementInfo::output(XMLHandler& xmlout) const {
  xmlout.set_root("KElementInfo");
  xmlout.put_child("JTimesTwo", make_string(getJtimestwo()));
  XMLHandler xmli;
  getRow().output(xmli);
  xmlout.put_child(xmli);
  getColumn().output(xmli);
  xmlout.put_child(xmli);
}

// ****************************************************************

KFitParamInfo::KFitParamInfo(XMLHandler& xmlin) {
  XMLHandler xmlk(xmlin, "KFitParamInfo");
  uint count1 = xmlk.count_among_children("PolynomialTerm");
  uint count2 = xmlk.count_among_children("PoleEnergy");
  uint count3 = xmlk.count_among_children("PoleCoupling");
  uint count4 = xmlk.count_among_children("StringExpressionParameter");
  if ((count1 + count2 + count3 + count4) != 1)
    throw(std::invalid_argument("Invalid KFitParamInfo XML"));
  if (count1 == 1) {
    XMLHandler xmlt(xmlk, "PolynomialTerm");
    uint power;
    xmlreadchild(xmlt, "Power", power);
    KElementInfo keleminfo(xmlt);
    set_poly_power(keleminfo, power);
  } else if (count2 == 1) {
    XMLHandler xmlt(xmlk, "PoleEnergy");
    uint poleindex, Jtimestwo;
    xmlreadchild(xmlt, "Index", poleindex);
    xmlreadchild(xmlt, "JTimesTwo", Jtimestwo);
    set_pole_energy(poleindex, Jtimestwo);
  } else if (count3 == 1) {
    XMLHandler xmlt(xmlk, "PoleCoupling");
    uint poleindex, Jtimestwo;
    xmlreadchild(xmlt, "Index", poleindex);
    xmlreadchild(xmlt, "JTimesTwo", Jtimestwo);
    KIndex kindex(xmlt);
    set_pole_coupling(kindex, poleindex, Jtimestwo);
  } else {
    XMLHandler xmlt(xmlk, "StringExpressionParameter");
    std::string param_name;
    xmlreadchild(xmlt, "ParameterName", param_name);
    KElementInfo keleminfo(xmlt);
    // Create a hash of the parameter name
    uint param_hash = 0;
    for (char c : param_name) {
      param_hash = param_hash * 31 + static_cast<uint>(c);
    }
    set_string_expr_param(keleminfo, param_hash);
  }
}

//    m_store1: left-most 4 bits encode the type of parameter
//         0 = term in a polynomial
//         1 = pole energy
//         2 = coupling in pole residue

//    if polynomial:
//         - rightmost 28 bits of m_store1 contains "n"
//               where parameter is coefficient of Ecm^n
//         - m_store2 contains the m_store of the K-matrix element
//               the polynomial is associated with

void KFitParamInfo::set_poly_power(const KElementInfo& keleminfo, uint power) {
  if (power > 32)
    throw(std::invalid_argument("Unsupported power in KFitParamInfo"));
  m_store1 = power;
  m_store2 = keleminfo.m_store;
}

//    if pole energy:
//         - rightmost 28 bits of m_store1 contain the integer index
//           identifying the pole (14 bits) and total J times two (14 bits)
//         - m_store2 = 0

void KFitParamInfo::set_pole_energy(uint pole_index, uint Jtimestwo) {
  if ((pole_index > 32) || (Jtimestwo > 128))
    throw(
        std::invalid_argument("Unsupported pole information in KFitParamInfo"));
  m_store1 = 1;
  m_store1 <<= 14;
  m_store1 |= pole_index;
  m_store1 <<= 14;
  m_store1 |= Jtimestwo;
  m_store2 = 0;
}

//    if coupling in pole residue
//         - rightmost 28 bits of m_store1 contain integer index
//              identify the pole (14 bits) and total J times two (14 bits)
//         - m_store2 contains m_store of the KIndex associated
//              with this parameter

void KFitParamInfo::set_pole_coupling(const KIndex& kindex, uint pole_index,
                                      uint Jtimestwo) {
  if ((pole_index > 32) || (Jtimestwo > 128))
    throw(
        std::invalid_argument("Unsupported pole information in KFitParamInfo"));
  m_store1 = 2;
  m_store1 <<= 14;
  m_store1 |= pole_index;
  m_store1 <<= 14;
  m_store1 |= Jtimestwo;
  m_store2 = kindex.m_store;
}

//    if string expression parameter:
//         - rightmost 28 bits of m_store1 contain a hash of the parameter name
//         - m_store2 contains the m_store of the K-matrix element
//               the parameter is associated with

void KFitParamInfo::set_string_expr_param(const KElementInfo& keleminfo, uint param_hash) {
  m_store1 = 3;
  m_store1 <<= 28;
  m_store1 |= (param_hash & 0xFFFFFFFu);
  m_store2 = keleminfo.m_store;
}

string KFitParamInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string KFitParamInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void KFitParamInfo::output(XMLHandler& xmlout) const {
  xmlout.set_root("KFitParamInfo");
  uint type = m_store1 >> 28;
  if (type == 0) {
    uint power = m_store1 & 0xFFFFFFFu;
    KElementInfo keleminfo(m_store2);
    XMLHandler xmlt("PolynomialTerm");
    xmlt.put_child("Power", make_string(power));
    XMLHandler xmlk;
    keleminfo.output(xmlk);
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);
  } else if (type == 1) {
    uint Jtimestwo = m_store1 & 0x3FFFu;
    uint poleindex = (m_store1 >> 14) & 0x3FFFu;
    XMLHandler xmlt("PoleEnergy");
    xmlt.put_child("Index", make_string(poleindex));
    xmlt.put_child("JTimesTwo", make_string(Jtimestwo));
    xmlout.put_child(xmlt);
  } else if (type == 2) {
    uint Jtimestwo = m_store1 & 0x3FFFu;
    uint poleindex = (m_store1 >> 14) & 0x3FFFu;
    XMLHandler xmlt("PoleCoupling");
    xmlt.put_child("Index", make_string(poleindex));
    xmlt.put_child("JTimesTwo", make_string(Jtimestwo));
    KIndex kindex(m_store2);
    XMLHandler xmlk;
    kindex.output(xmlk);
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);
  } else if (type == 3) {
    uint param_hash = m_store1 & 0xFFFFFFFu;
    KElementInfo keleminfo(m_store2);
    XMLHandler xmlt("StringExpressionParameter");
    xmlt.put_child("ParameterHash", make_string(param_hash));
    XMLHandler xmlk;
    keleminfo.output(xmlk);
    xmlt.put_child(xmlk);
    xmlout.put_child(xmlt);
  }
}

string KFitParamInfo::getMCObsName() const {
  uint type = m_store1 >> 28;
  if (type == 0) {
    uint power = m_store1 & 0xFFFFFFFu;
    KElementInfo keleminfo(m_store2);
    return string("KPoly") + make_stringtuple(power, keleminfo.getJtimestwo()) +
           make_stringtuple(keleminfo.getRowL(), keleminfo.getRowStimestwo(),
                            keleminfo.getRowChannelIndex()) +
           make_stringtuple(keleminfo.getColumnL(),
                            keleminfo.getColumnStimestwo(),
                            keleminfo.getColumnChannelIndex());
  } else if (type == 1) {
    uint Jtimestwo = m_store1 & 0x3FFFu;
    uint poleindex = (m_store1 >> 14) & 0x3FFFu;
    return string("KPoleE") + make_stringtuple(poleindex, Jtimestwo);
  } else if (type == 2) {
    uint Jtimestwo = m_store1 & 0x3FFFu;
    uint poleindex = (m_store1 >> 14) & 0x3FFFu;
    KIndex kindex(m_store2);
    return string("KPoleC") + make_stringtuple(poleindex, Jtimestwo) +
           make_stringtuple(kindex.getL(), kindex.getStimestwo(),
                            kindex.getChannelIndex());
  } else if (type == 3) {
    // String expression parameter
    uint param_hash = m_store1 & 0xFFFFFFFu;
    KElementInfo keleminfo(m_store2);
    return string("KStrExpr") + make_stringtuple(param_hash, keleminfo.getJtimestwo()) +
           make_stringtuple(keleminfo.getRowL(), keleminfo.getRowStimestwo(),
                            keleminfo.getRowChannelIndex()) +
           make_stringtuple(keleminfo.getColumnL(),
                            keleminfo.getColumnStimestwo(),
                            keleminfo.getColumnChannelIndex());
  }
  return string("");
}

void KFitParamInfo::setFromMCObsName(const std::string& obsname) {
  try {
    if (obsname.substr(0, 5) == string("KPoly")) {
      vector<string> tuples(extract_stringtuples(obsname));
      if (tuples.size() != 3)
        throw(std::invalid_argument("Error"));
      int power, Jtimestwo, L, Stimestwo, a, Lp, Stimestwop, ap;
      read_stringtuple(tuples[0], power, Jtimestwo);
      read_stringtuple(tuples[1], L, Stimestwo, a);
      read_stringtuple(tuples[2], Lp, Stimestwop, ap);
      KElementInfo keleminfo(Jtimestwo, L, Stimestwo, a, Lp, Stimestwop, ap);
      set_poly_power(keleminfo, power);
    } else if (obsname.substr(0, 6) == string("KPoleE")) {
      vector<string> tuples(extract_stringtuples(obsname));
      if (tuples.size() != 1)
        throw(std::invalid_argument("Error"));
      int poleindex, Jtimestwo;
      read_stringtuple(tuples[0], poleindex, Jtimestwo);
      set_pole_energy(poleindex, Jtimestwo);
    } else if (obsname.substr(0, 6) == string("KPoleC")) {
      vector<string> tuples(extract_stringtuples(obsname));
      if (tuples.size() != 2)
        throw(std::invalid_argument("Error"));
      int poleindex, Jtimestwo, L, Stimestwo, a;
      read_stringtuple(tuples[0], poleindex, Jtimestwo);
      read_stringtuple(tuples[1], L, Stimestwo, a);
      KIndex kindex(L, Stimestwo, a);
      set_pole_coupling(kindex, poleindex, Jtimestwo);
    } else if (obsname.substr(0, 8) == string("KStrExpr")) {
      vector<string> tuples(extract_stringtuples(obsname));
      if (tuples.size() != 3)
        throw(std::invalid_argument("Error"));
      int param_hash, Jtimestwo, L, Stimestwo, a, Lp, Stimestwop, ap;
      read_stringtuple(tuples[0], param_hash, Jtimestwo);
      read_stringtuple(tuples[1], L, Stimestwo, a);
      read_stringtuple(tuples[2], Lp, Stimestwop, ap);
      KElementInfo keleminfo(Jtimestwo, L, Stimestwo, a, Lp, Stimestwop, ap);
      set_string_expr_param(keleminfo, param_hash);
    } else
      throw(std::invalid_argument("Error"));
  } catch (const std::exception& xp) {
    throw(std::invalid_argument("Could not setFromMCObsName in KFitParamInfo"));
  }
}

// ****************************************************************

DecayChannelInfo::DecayChannelInfo(const std::string& particle1_name,
                                   const std::string& particle2_name,
                                   uint spin1_times_two, uint spin2_times_two,
                                   bool identical, bool same_intrinsic_parities)
    : m_name1(tidyString(particle1_name)), m_name2(tidyString(particle2_name)) {
  check_encode(spin1_times_two, spin2_times_two, identical,
               same_intrinsic_parities);
}

void DecayChannelInfo::check_encode(uint spin1_times_two, uint spin2_times_two,
                                    bool identical,
                                    bool same_intrinsic_parities) {
  if ((spin1_times_two >= 32768) || (spin2_times_two >= 32768) ||
      (m_name1.empty()) || (m_name2.empty()))
    throw(std::invalid_argument(
        "Invalid DecayChannelInfo: empty name, or spin too large"));
  if ((identical) &&
      ((m_name1 != m_name2) || (spin1_times_two != spin2_times_two)))
    throw(std::invalid_argument(
        "DecayChannelInfo identical particles but different properties"));
  m_store = spin1_times_two;
  m_store <<= 15;
  m_store |= spin2_times_two;
  m_store <<= 1;
  if (identical)
    m_store |= 1u;
  m_store <<= 1;
  if (same_intrinsic_parities)
    m_store |= 1u;
}

DecayChannelInfo::DecayChannelInfo(XMLHandler& xmlin) {
  uint spin1_times_two, spin2_times_two;
  bool identical, same_intrinsic_parities;
  try {
    XMLHandler xmlinfo(xmlin, "DecayChannelInfo");
    xmlreadchild(xmlinfo, "Particle1Name", m_name1);
    xmlreadchild(xmlinfo, "Spin1TimesTwo", spin1_times_two);
    if (xml_tag_count(xmlinfo, "Identical") > 0) {
      string name2(m_name1);
      xmlreadifchild(xmlinfo, "Particle2Name", name2);
      uint s2 = spin1_times_two;
      xmlreadifchild(xmlinfo, "Spin2TimesTwo", s2);
      string ipprod("same");
      xmlreadifchild(xmlinfo, "IntrinsicParities", ipprod);
      if ((name2 != m_name1) || (s2 != spin1_times_two) || (ipprod != "same"))
        throw(std::invalid_argument(
            string("Input not consistent with identical particles")));
      m_name2 = m_name1;
      spin2_times_two = spin1_times_two;
      identical = true;
      same_intrinsic_parities = true;
    } else {
      identical = false;
      xmlreadchild(xmlinfo, "Particle2Name", m_name2);
      xmlreadchild(xmlinfo, "Spin2TimesTwo", spin2_times_two);
      string reply;
      xmlreadchild(xmlinfo, "IntrinsicParities", reply);
      if (reply == "same")
        same_intrinsic_parities = true;
      else if (reply == "opposite")
        same_intrinsic_parities = false;
      else
        throw(std::invalid_argument("Invalid IntrinsicParities tag"));
    }
    check_encode(spin1_times_two, spin2_times_two, identical,
                 same_intrinsic_parities);
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Invalid input to DecayChannelInfo: ") +
                                xp.what()));
  }
}

DecayChannelInfo::DecayChannelInfo(const DecayChannelInfo& dcinfo)
    : m_name1(dcinfo.m_name1), m_name2(dcinfo.m_name2),
      m_store(dcinfo.m_store) {}

DecayChannelInfo& DecayChannelInfo::operator=(const DecayChannelInfo& dcinfo) {
  m_name1 = dcinfo.m_name1;
  m_name2 = dcinfo.m_name2;
  m_store = dcinfo.m_store;
  return *this;
}

string DecayChannelInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string DecayChannelInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void DecayChannelInfo::output(XMLHandler& xmlout) const {
  xmlout.set_root("DecayChannelInfo");
  xmlout.put_child("Particle1Name", m_name1);
  xmlout.put_child("Spin1TimesTwo", make_string(getSpin1timestwo()));
  if (areIdentical())
    xmlout.put_child("Identical");
  else {
    xmlout.put_child("Particle2Name", m_name2);
    xmlout.put_child("Spin2TimesTwo", make_string(getSpin2timestwo()));
    if (sameIntrinsicParities())
      xmlout.put_child("IntrinsicParities", "same");
    else
      xmlout.put_child("IntrinsicParities", "opposite");
  }
}

ostream& operator<<(ostream& os, const DecayChannelInfo& dcinfo) {
  os << dcinfo.output();
  return os;
}

// *******************************************************************
