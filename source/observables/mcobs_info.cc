#include "mcobs_info.h"
#include "encoder.h"
#include "multi_compare.h"

using namespace std;

// ***************************************************************

MCObsInfo::MCObsInfo() : icode(1) {
  icode[0] = 0; // default is zero particles
}

MCObsInfo::MCObsInfo(XMLHandler& xml_in) {
  try {
    set<string> tags;
    tags.insert("MCObservable");
    tags.insert("MCObs");
    ArgsHandler xin(xml_in, tags);
    string rtag = xin.getInputRootTag();
    if (rtag == "MCObservable") {
      assign(xin);
      return;
    } else if (rtag == "MCObs") {
      assign_from_string(xin.getString("MCObs"));
      return;
    } else
      throw(std::invalid_argument("Invalid XML input"));
  } catch (const std::exception& errmsg) {
    throw(std::invalid_argument(string("MCObsInfo construction failed: \n") +
                                string(errmsg.what()) + string("\nInput XML:") +
                                xml_in.output()));
  }
}

void MCObsInfo::assign(ArgsHandler& xt) {
  string name(xt.getString("ObsName"));
  uint index = 0;
  xt.getOptionalUInt("Index", index);
  bool simple = xt.queryTag("Simple");
  ComplexArg arg = RealPart;
  string reply = "Re";
  xt.getOptionalString("Arg", reply);
  if ((reply == "ImaginaryPart") || (reply == "Im"))
    arg = ImaginaryPart;
  else if ((reply == "RealPart") || (reply == "Re"))
    arg = RealPart;
  else
    throw(std::invalid_argument("Invalid Arg tag"));
  encode(name, index, simple, arg);
}

void MCObsInfo::assign_from_string(const string& opstring) {
  string opstr(tidyString(opstring));
  vector<string> tokens = split(opstr, ' ');
  if ((tokens.size() > 4) && (tokens.size() < 1))
    throw(std::runtime_error(""));
  string name(tokens[0]);
  uint index = 0;
  if (tokens.size() > 1) {
    extract_from_string(tokens[1], index);
  }
  bool simple;
  if (tokens.size() < 3)
    simple = false;
  else if (tokens[2] == "s")
    simple = true;
  else if (tokens[2] == "n")
    simple = false;
  else
    throw(std::runtime_error(""));
  ComplexArg arg;
  if (tokens.size() < 4)
    arg = RealPart;
  else if (tokens[3] == "re")
    arg = RealPart;
  else if (tokens[3] == "im")
    arg = ImaginaryPart;
  else
    throw(std::runtime_error(""));
  encode(name, index, simple, arg);
}

vector<string> MCObsInfo::split(const string& astr, char delimiter) const {
  vector<string> tokens;
  size_t lastpos = astr.find_first_not_of(delimiter);
  size_t pos = (lastpos == string::npos)
                   ? string::npos
                   : astr.find_first_of(delimiter, lastpos + 1);
  while (lastpos != string::npos) {
    if (pos == string::npos)
      pos = astr.length();
    tokens.push_back(astr.substr(lastpos, pos - lastpos));
    lastpos = astr.find_first_not_of(delimiter, pos + 1);
    pos = (lastpos == string::npos)
              ? string::npos
              : astr.find_first_of(delimiter, lastpos + 1);
  }
  return tokens;
}

MCObsInfo::MCObsInfo(const string& obsname, uint index, bool simple,
                     ComplexArg arg) {
  encode(obsname, index, simple, arg);
}

void MCObsInfo::setToRealPart() { set_real_part(); }

void MCObsInfo::setToImaginaryPart() {
  if (isVacuum())
    throw(std::invalid_argument("Cannot set vacuum to imaginary part"));
  set_imag_part();
}

void MCObsInfo::resetObsIndex(uint ind) {
  if (isPrimary())
    throw(std::invalid_argument("Cannot reset index for primary observable"));
  set_index(ind);
}

bool MCObsInfo::isVacuum() const { return (icode[0] == 0u); }

bool MCObsInfo::isRealPart() const { return ((icode[0] & 1u) == 0); }

bool MCObsInfo::isImaginaryPart() const { return ((icode[0] & 1u) == 1); }

bool MCObsInfo::isSimple() const { return (((icode[0] >> 1) & 1u) == 0); }

bool MCObsInfo::isNonSimple() const { return (((icode[0] >> 1) & 1u) != 0); }

bool MCObsInfo::isPrimary() const { return (((icode[0] >> 2) & 1u) == 0); }

bool MCObsInfo::isSecondary() const { return (((icode[0] >> 2) & 1u) != 0); }

string MCObsInfo::getObsName() const {
  if (isPrimary()) {
    throw(std::invalid_argument("getObsName called for primary observable"));
  }
  return get_obs_name();
}

uint MCObsInfo::getObsIndex() const {
  if (isPrimary()) {
    throw(std::invalid_argument("getObsName called for primary observable"));
  }
  return get_obs_index();
}

string MCObsInfo::output(bool longform, int indent) const {
  XMLHandler xmlout;
  output(xmlout, longform);
  return xmlout.output(indent);
}

string MCObsInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void MCObsInfo::output(XMLHandler& xmlout, bool longform) const {
  xmlout.set_root("MCObservable");
  if (isVacuum()) {
    xmlout.put_child("Vacuum");
  } else {
    string obsname = get_obs_name();
    uint index = get_obs_index();
    if (longform) {
      xmlout.put_child("ObsName", obsname);
      xmlout.put_child("Index", make_string(index));
      if (isSimple())
        xmlout.put_child("Simple");
      if (isRealPart())
        xmlout.put_child("Arg", "RealPart");
      else
        xmlout.put_child("Arg", "ImaginaryPart");
    } else {
      string infostr(obsname);
      infostr += " " + make_string(index);
      if (isSimple())
        infostr += " s";
      else
        infostr += " n";
      infostr += (isRealPart() ? " re" : " im");
      xmlout.set_root("MCObs", infostr);
    }
  }
}

bool MCObsInfo::operator==(const MCObsInfo& rhs) const {
  return multiEqual(icode, rhs.icode);
}

bool MCObsInfo::operator!=(const MCObsInfo& rhs) const {
  return multiNotEqual(icode, rhs.icode);
}

bool MCObsInfo::operator<(const MCObsInfo& rhs) const {
  return multiLessThan(icode, rhs.icode);
}

//  private routines

void MCObsInfo::encode(const string& obsname, uint index, bool simple,
                       ComplexArg arg) {
  uint nchar = obsname.length();
  if (nchar > 64) {
    throw(std::invalid_argument(
        "MCObsInfo name cannot be longer than 64 characters"));
  }
  vector<uint> namecode;
  encode_string_to_uints(obsname, 64, namecode);
  icode.resize(namecode.size() + 1);
  std::copy(namecode.begin(), namecode.end(), icode.begin() + 1);
  uint tcode = index;
  tcode <<= 1;
  tcode |= 1u; // secondary
  tcode <<= 1;
  if (!simple)
    tcode |= 1u;
  tcode <<= 1;
  if (arg == ImaginaryPart)
    tcode |= 1u;
  icode[0] = tcode;
}

void MCObsInfo::set_real_part() {
  icode[0] &= ~1u; // clear the bit
}

void MCObsInfo::set_imag_part() {
  icode[0] |= 1u; // set the bit
}

void MCObsInfo::set_index(uint index) {
  uint tcode = index;
  tcode <<= 3;
  tcode |= (icode[0] & 7u);
  icode[0] = tcode;
}

void MCObsInfo::set_arg(ComplexArg arg) {
  if (arg == RealPart)
    set_real_part();
  else
    set_imag_part();
}

string MCObsInfo::get_obs_name() const {
  vector<uint> namecode(icode.begin() + 1, icode.end());
  return decode_uints_to_string(namecode);
}

uint MCObsInfo::get_obs_index() const { return (icode[0] >> 3); }

void MCObsInfo::copyTo(unsigned int* buf) const {
  buf[0] = icode.size();
  if (buf[0] >= max_ints)
    throw(std::runtime_error(
        "icode in MCObsInfo too large for record key output"));
  for (uint i = 0; i < buf[0]; i++)
    buf[i + 1] = icode[i];
  for (uint i = buf[0] + 1; i < max_ints; i++)
    buf[i] = 0;
}

MCObsInfo::MCObsInfo(const unsigned int* buf) {
  if (buf[0] >= max_ints)
    throw(std::runtime_error(
        "icode in MCObsInfo too large for record key output"));
  icode.resize(buf[0]);
  std::copy(buf + 1, buf + 1 + buf[0], icode.begin());
}

// ******************************************************************************
