#include "sampling_info.h"
#include "args_handler.h"
#include <stdexcept>
using namespace std;

// *************************************************************

MCSamplingInfo::MCSamplingInfo() : icode(0), icode2(100) {}

MCSamplingInfo::MCSamplingInfo(XMLHandler& xml_in) {
  try {
    ArgsHandler xmlr(xml_in, "MCSamplingInfo");
    bool jack = xmlr.queryTag("Jackknife");
    bool boot = xmlr.queryTag("Bootstrapper");
    if (jack == boot) {
      throw(std::invalid_argument("Invalid mode in MCSamplingInfo"));
    }
    if (jack) {
      icode = 0;
      uint nbins = 100;
      ArgsHandler xmlj(xmlr, "Jackknife");
      xmlj.getOptionalUInt("NumberBins", nbins);
      icode2 = nbins;
      return;
    }
    ArgsHandler xmlb(xmlr, "Bootstrapper");
    uint num_resamplings = 1024;
    uint bootseed = 0;
    uint bootskip = 64;
    xmlb.getOptionalUInt("NumberResamplings", num_resamplings);
    xmlb.getOptionalUInt("Seed", bootseed);
    xmlb.getOptionalUInt("BootSkip", bootskip);
    encode_bootstrap(num_resamplings, bootseed, bootskip);
  } catch (std::exception& xp) {
    throw(std::invalid_argument(string("MCSamplingInfo creation failed: ") +
                                xp.what()));
  }
}

MCSamplingInfo::MCSamplingInfo(XMLHandler& xml_in, const MCBinsInfo& bin_info) {
  try {
    ArgsHandler xmlr(xml_in, "MCSamplingInfo");
    bool jack = xmlr.queryTag("Jackknife");
    bool boot = xmlr.queryTag("Bootstrapper");
    if (jack == boot) {
      throw(std::invalid_argument("Invalid mode in MCSamplingInfo"));
    }
    if (jack) {
      uint njackbins = bin_info.getNumberOfBins();
      icode = 0;
      icode2 = njackbins;
      uint nbinsread = 0;
      ArgsHandler xmlj(xmlr, "Jackknife");
      xmlj.getOptionalUInt("NumberBins", nbinsread);
      if ((nbinsread != 0) && (nbinsread != njackbins))
        throw(
            std::invalid_argument("Invalid number of bins in MCSamplingInfo"));
      return;
    }
    ArgsHandler xmlb(xmlr, "Bootstrapper");
    uint num_resamplings = 1024;
    uint bootseed = 0;
    uint bootskip = 64;
    xmlb.getOptionalUInt("NumberResamplings", num_resamplings);
    xmlb.getOptionalUInt("Seed", bootseed);
    xmlb.getOptionalUInt("BootSkip", bootskip);
    encode_bootstrap(num_resamplings, bootseed, bootskip);
  } catch (std::exception& xp) {
    throw(std::invalid_argument(string("MCSamplingInfo creation failed: ") +
                                xp.what()));
  }
}

MCSamplingInfo::MCSamplingInfo(const MCSamplingInfo& fin)
    : icode(fin.icode), icode2(fin.icode2) {}

MCSamplingInfo& MCSamplingInfo::operator=(const MCSamplingInfo& fin) {
  icode = fin.icode;
  icode2 = fin.icode2;
  return *this;
}

MCSamplingInfo::MCSamplingInfo(uint nbootsamp, unsigned long bootseed,
                               uint bootskip) {
  encode_bootstrap(nbootsamp, bootseed, bootskip);
}

MCSamplingInfo::MCSamplingInfo(uint njackbins) {
  icode = 0;
  icode2 = njackbins;
}

bool MCSamplingInfo::isJackknifeMode() const { return (icode == 0); }

bool MCSamplingInfo::isBootstrapMode() const { return (icode != 0); }

SamplingMode MCSamplingInfo::getSamplingMode() const {
  return (icode == 0) ? Jackknife : Bootstrap;
}

unsigned int MCSamplingInfo::getNumberOfReSamplings() const {
  return (icode == 0) ? icode2 : icode & 0xFFFFFu;
}

// private member; used by KBObsHandler

unsigned int MCSamplingInfo::getNumberOfBootstrapReSamplings() const {
  return icode & 0xFFFFFu;
}

unsigned long
MCSamplingInfo::getRNGSeed() const // returns zero for jackknife mode
{
  return (icode != 0) ? icode2 : 0;
}

unsigned int
MCSamplingInfo::getSkipValue() const // returns zero for jackknife mode
{
  return (icode >> 20);
}

void MCSamplingInfo::checkEqual(const MCSamplingInfo& rhs) const {
  if ((icode != rhs.icode) || (icode2 != rhs.icode2)) {
    cerr << "MCSamplingInfo checkEqual failed" << endl;
    cerr << "LHS:" << endl
         << output() << endl
         << "RHS:" << endl
         << rhs.output() << endl;
    throw(std::invalid_argument("MCSamplingInfo checkEqual failed"));
  }
}

bool MCSamplingInfo::operator==(const MCSamplingInfo& rhs) const {
  return ((icode == rhs.icode) && (icode2 == rhs.icode2));
}

bool MCSamplingInfo::operator!=(const MCSamplingInfo& rhs) const {
  return ((icode != rhs.icode) || (icode2 != rhs.icode2));
}

void MCSamplingInfo::output(XMLHandler& xmlout) const {
  xmlout.set_root("MCSamplingInfo");
  if (icode == 0) {
    XMLHandler xmlj("Jackknife");
    xmlj.put_child("NumberBins", make_string(icode2));
    xmlout.put_child(xmlj);
    return;
  }
  XMLHandler xmlb("Bootstrapper");
  uint nsamp = icode & 0xFFFFFu;
  uint skip = icode >> 20;
  xmlb.put_child("NumberResamplings", make_string(nsamp));
  xmlb.put_child("Seed", make_string(icode2));
  xmlb.put_child("BootSkip", make_string(skip));
  xmlout.put_child(xmlb);
}

string MCSamplingInfo::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

string MCSamplingInfo::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

void MCSamplingInfo::encode_bootstrap(unsigned int num_resamplings,
                                      unsigned long seed,
                                      unsigned int skip_value) {
  if (num_resamplings < 1)
    throw(std::invalid_argument(
        "Number of resamplings > 0 required in MCSamplingInfo"));
  if (num_resamplings >= 1048576)
    throw(std::invalid_argument(
        "Number of resamplings too large in MCSamplingInfo"));
  if (skip_value >= 4096)
    throw(std::invalid_argument("Skip value too large in MCSamplingInfo"));
  icode = skip_value;
  icode <<= 20;
  icode |= num_resamplings;
  icode2 = seed;
}

// ***************************************************************
