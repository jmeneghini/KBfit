#include "obs_get_handler.h"
#include "args_handler.h"
#include "io_handler_fstream.h"
#include "io_handler_hdf5.h"

using namespace std;

// *************************************************************************

MCObsGetHandler::MCObsGetHandler(XMLHandler& xmlin) {
  m_sampinfo = 0;
  XMLHandler xmlr(xmlin, "KBObservables");
  xml_tag_assert(xmlr, "KBObservables");
  try {
    MCSamplingInfo sampinfo(xmlr);
    initialize(sampinfo);
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(
        "Could not read MCSamplingInfo in MCObsGetHandler"));
  }

  bool verbose = false;
  if (xml_tag_count(xmlr, "Verbose") > 0)
    verbose = true;

  //  read the input file lists to find the data

  set<string> sampfiles;
  if (xmlr.query_unique_to_among_children("SamplingData")) {
    XMLHandler xmlp(xmlr, "SamplingData");
    {
      list<XMLHandler> infxml;
      infxml = xmlp.find("FileName");
      for (list<XMLHandler>::iterator it = infxml.begin(); it != infxml.end();
           ++it) {
        ArgsHandler xmls(*it);
        string fname(xmls.getString("FileName"));
        sampfiles.insert(fname);
      }
    }
  }
  if (sampfiles.empty()) {
    m_logger << "No sampling files given" << endl;
    return;
  }

  //  read requested observables (optionally given)

  std::map<MCEnsembleInfo, std::set<MCObsInfo>> obs_map;
  int speccount = xml_tag_count(xmlr, "Specifications");
  if (speccount > 1) {
    throw(std::invalid_argument("Multiple <Specifications> tags not allowed"));
  }
  try {
    XMLHandler xmlm(xmlin, "Specifications");
    list<XMLHandler> ensxml;
    ensxml = xmlm.find_among_children("ObsSamplings");
    for (list<XMLHandler>::iterator it = ensxml.begin(); it != ensxml.end();
         ++it) {
      MCEnsembleInfo mcens(*it);
      if (obs_map.find(mcens) != obs_map.end())
        throw(std::runtime_error("Multiple instances of same ensemble tag"));
      if (xml_tag_count(*it, "All") > 0) {
        set<MCObsInfo> mcobsset;
        obs_map.insert(make_pair(mcens, mcobsset));
      } else {
        list<XMLHandler> obsxml;
        set<MCObsInfo> mcobsset;
        obsxml = it->find_among_children("MCObservable");
        for (list<XMLHandler>::iterator mt = obsxml.begin(); mt != obsxml.end();
             ++mt) {
          MCObsInfo obskey(*mt);
          mcobsset.insert(obskey);
        }
        if (mcobsset.size() > 0)
          obs_map.insert(make_pair(mcens, mcobsset));
      }
    }
  } catch (const std::exception& xp) {
  }

  connectSamplingFiles(sampfiles, obs_map, verbose);
}

void MCObsGetHandler::initialize(const MCSamplingInfo& sampinfo) {
  if ((m_sampinfo != 0) && ((*m_sampinfo) == sampinfo))
    return;
  m_sampinfo = new MCSamplingInfo(sampinfo);
  m_sampmode = m_sampinfo->getSamplingMode();
  //  add "indep" to m_bininfos
  MCEnsembleInfo mcindep(m_sampinfo->getNumberOfReSamplings());
  m_bininfos.insert(make_pair(mcindep, MCBinsInfo(mcindep)));
}

MCObsGetHandler::MCObsGetHandler(const MCSamplingInfo& sampinfo) {
  m_sampinfo = 0;
  initialize(sampinfo);
}

MCObsGetHandler::MCObsGetHandler(const MCSamplingInfo& sampinfo,
                                 const std::set<std::string>& sampfiles,
                                 bool verbose) {
  m_sampinfo = 0;
  initialize(sampinfo);
  connectSamplingFiles(sampfiles, verbose);
}

MCObsGetHandler::MCObsGetHandler(
    const MCSamplingInfo& sampinfo, const std::set<std::string>& sampfiles,
    const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep,
    bool verbose) {
  m_sampinfo = 0;
  initialize(sampinfo);
  connectSamplingFiles(sampfiles, keys_to_keep, verbose);
}

MCObsGetHandler::~MCObsGetHandler() {
  clear();
  delete m_sampinfo;
}

void MCObsGetHandler::connectSamplingFiles(const std::set<string>& sampfiles,
                                           bool verbose) {
  std::map<MCEnsembleInfo, std::set<MCObsInfo>> keys_to_keep;
  connectSamplingFiles(sampfiles, keys_to_keep, verbose); // keeps all
}

void MCObsGetHandler::connectSamplingFiles(
    const std::set<std::string>& sampfiles,
    const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep,
    bool verbose) {
  if (verbose)
    m_logger << endl << "Observable get handler connecting to files:" << endl;
  // now connect the files
  for (set<string>::iterator st = sampfiles.begin(); st != sampfiles.end();
       st++) {
    if (verbose)
      m_logger << "Now connecting file " << *st << endl;
    connect_a_samplings_file(*st, keys_to_keep);
  }

  // now check that all requested observables are available
  if (!(keys_to_keep.empty())) {
    try {
      for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator pt =
               m_sampsdh.begin();
           pt != m_sampsdh.end(); ++pt) {
        const MCEnsembleInfo& mcens = pt->first;
        std::map<MCEnsembleInfo, set<MCObsInfo>>::const_iterator kt =
            keys_to_keep.find(mcens);
        if ((kt != keys_to_keep.end()) && (!((kt->second).empty()))) {
          if (!(pt->second->keepKeys(kt->second)))
            throw(std::runtime_error(" The observable " +
                                     kt->second.begin()->output() +
                                     " from the ensemble " + mcens.output()));
        }
      }
    } catch (const std::exception& xp) {
      string errmsg("Error");
      errmsg += xp.what();
      errmsg += getCurrentLog().str();
      clear();
      throw(std::runtime_error(errmsg));
    }
  }

  if (verbose) {
    XMLHandler xmlfmap;
    getFileMap(xmlfmap);
    m_logger << endl << "File map:" << xmlfmap.output() << endl;
    for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator pt =
             m_sampsdh.begin();
         pt != m_sampsdh.end(); ++pt) {
      const MCEnsembleInfo& mcens = pt->first;
      m_logger << "MCEnsemble connection:" << mcens.output()
               << " Size = " << pt->second->size() << endl;
      if (verbose) {
        XMLHandler xmlkeys;
        pt->second->outputKeys(xmlkeys);
        m_logger << "Keys:" << endl << xmlkeys.output() << endl;
      }
    }
  }
}

// if map "keys_to_keep" is empty:
//    -- keep all keys encountered for ALL ensembles encountered
// else if "keys_to_keep" DOES have an entry for an ensemble:
//    -- if entry is empty set, keep all keys for that ensemble
//    -- if entry is nonempty set, keep only those keys in the set
// else if "keys_to_keep" does NOT have an entry for an ensemble:
//    -- no keys are to be kept for that ensemble

void MCObsGetHandler::connect_a_samplings_file(
    const std::string& file_name,
    const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep) {
  try {
    std::pair<MCBinsInfo, MCSamplingInfo> finfo = get_info_from_file(file_name);
    if (finfo.second != (*m_sampinfo)) {
      m_logger << "Data in " << file_name
               << " does not match requested sampling info; skipping" << endl;
      return;
    }
    MCEnsembleInfo mcens = finfo.first.getMCEnsembleInfo();
    if ((!(keys_to_keep.empty())) &&
        (keys_to_keep.find(mcens) == keys_to_keep.end()))
      return; // do not keep
    std::map<MCEnsembleInfo, MCBinsInfo>::iterator bt = m_bininfos.find(mcens);
    if (bt == m_bininfos.end()) {
      m_bininfos.insert(make_pair(mcens, finfo.first));
    } else if ((bt->second) != finfo.first) {
      m_logger
          << "Bin info must match for all files associated with same ensemble"
          << endl;
      throw(std::runtime_error(
          "Bin info must match for all files associated with same ensemble"));
    }
    std::map<MCEnsembleInfo, SamplingsGetHandler*>::iterator it =
        m_sampsdh.find(mcens);
    if (it == m_sampsdh.end()) {
      set<string> sampfiles;
      SamplingsGetHandler* sgptr =
          new SamplingsGetHandler(finfo.first, finfo.second, sampfiles);
      pair<std::map<MCEnsembleInfo, SamplingsGetHandler*>::iterator, bool> ret =
          m_sampsdh.insert(make_pair(mcens, sgptr));
      if (!(ret.second))
        throw(std::runtime_error("Error inserting into sampling handler map"));
      it = ret.first;
    }
    std::map<MCEnsembleInfo, std::set<MCObsInfo>>::const_iterator kt =
        keys_to_keep.find(mcens);
    std::set<MCObsInfo> dummy; // an empty set
    const std::set<MCObsInfo>& keepset =
        (kt != keys_to_keep.end()) ? kt->second : dummy;
    it->second->addFile(file_name, keepset);
  } catch (const std::exception& msg) {
    string errmsg("Error");
    errmsg += msg.what();
    errmsg += getCurrentLog().str();
    clear();
    throw(std::runtime_error(errmsg));
  }
}

void MCObsGetHandler::disconnectSamplingFiles(const set<string>& filenames) {
  for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator it =
           m_sampsdh.begin();
       it != m_sampsdh.end(); it++) {
    for (set<string>::const_iterator ft = filenames.begin();
         ft != filenames.end(); ++ft) {
      it->second->removeFile(*ft);
    }
  }
}

void MCObsGetHandler::disconnectAllSamplingFiles() {
  for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator it =
           m_sampsdh.begin();
       it != m_sampsdh.end(); it++)
    delete it->second;
  m_sampsdh.clear();
  m_bininfos.clear();
}

// clears everything except m_sampinfo, m_sampmode

void MCObsGetHandler::clear() {
  for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator it =
           m_sampsdh.begin();
       it != m_sampsdh.end(); it++)
    delete it->second;
  m_sampsdh.clear();
  m_bininfos.clear();
  clearLog();
}

// "setSamplingMode" does complete clear if new mode differs
// from current;  nothing happens if new=current mode

void MCObsGetHandler::setSamplingInfo(const MCSamplingInfo& sampinfo) {
  if ((*m_sampinfo) == sampinfo)
    return;
  clear();
  initialize(sampinfo);
}

const MCBinsInfo&
MCObsGetHandler::getBinsInfo(const MCEnsembleInfo& mcens) const {
  std::map<MCEnsembleInfo, MCBinsInfo>::const_iterator it =
      m_bininfos.find(mcens);
  if (it == m_bininfos.end())
    throw(std::runtime_error("Invalid MCEnsembleInfo in getBinsInfo"));
  return it->second;
}

const MCSamplingInfo& MCObsGetHandler::getSamplingInfo() const {
  if (m_sampinfo == 0)
    throw(std::runtime_error("getSamplingInfo failed"));
  return *m_sampinfo;
}

uint MCObsGetHandler::getNumberOfResamplings() const {
  return m_sampinfo->getNumberOfReSamplings();
}

SamplingMode MCObsGetHandler::getSamplingMode() const { return m_sampmode; }

bool MCObsGetHandler::isJackknifeMode() const {
  return (m_sampmode == Jackknife);
}

bool MCObsGetHandler::isBootstrapMode() const {
  return (m_sampmode == Bootstrap);
}

std::set<MCEnsembleInfo> MCObsGetHandler::getEnsembleInfos() const {
  set<MCEnsembleInfo> eset;
  for (map<MCEnsembleInfo, MCBinsInfo>::const_iterator it = m_bininfos.begin();
       it != m_bininfos.end(); ++it)
    eset.insert(it->first);
  return eset;
}

void MCObsGetHandler::getSamplings(const KBObsInfo& obsinfo,
                                   RVector& samplings) const {
  map<MCEnsembleInfo, SamplingsGetHandler*>::const_iterator it =
      m_sampsdh.find(obsinfo.getMCEnsembleInfo());
  if (it == m_sampsdh.end())
    throw(std::invalid_argument(
        string("getSamplings fails due to unavailable sampling for ") +
        obsinfo.str()));
  it->second->getData(obsinfo.getMCObsInfo(), samplings);
}

bool MCObsGetHandler::getSamplingsMaybe(const KBObsInfo& obsinfo,
                                        RVector& samplings) const {
  samplings.clear();
  map<MCEnsembleInfo, SamplingsGetHandler*>::const_iterator it =
      m_sampsdh.find(obsinfo.getMCEnsembleInfo());
  if (it == m_sampsdh.end()) {
    return false;
  }
  return it->second->getDataMaybe(obsinfo.getMCObsInfo(), samplings);
}

bool MCObsGetHandler::querySamplings(const KBObsInfo& obsinfo) const {
  map<MCEnsembleInfo, SamplingsGetHandler*>::const_iterator it =
      m_sampsdh.find(obsinfo.getMCEnsembleInfo());
  if (it == m_sampsdh.end())
    return false;
  return it->second->queryData(obsinfo.getMCObsInfo());
}

void MCObsGetHandler::close() {
  for (map<MCEnsembleInfo, SamplingsGetHandler*>::iterator it =
           m_sampsdh.begin();
       it != m_sampsdh.end(); it++)
    it->second->close();
}

void MCObsGetHandler::getFileMap(XMLHandler& xmlout) const {
  xmlout.set_root("FileMap");
  for (map<MCEnsembleInfo, SamplingsGetHandler*>::const_iterator it =
           m_sampsdh.begin();
       it != m_sampsdh.end(); it++) {
    XMLHandler xmls;
    it->first.output(xmls);
    XMLHandler xmle("EnsembleMapping");
    xmle.put_child(xmls);
    set<string> fnames(it->second->getFileNames());
    for (set<string>::iterator st = fnames.begin(); st != fnames.end(); st++) {
      xmle.put_child("FileName", *st);
    }
    xmlout.put_child(xmle);
  }
}

//   private members

std::pair<MCBinsInfo, MCSamplingInfo>
MCObsGetHandler::get_info_from_file(const std::string& filename) {
  // First, try to peek the file ID to determine the format and verify it's a
  // Sigmond file
  std::string ID;
  IOFSTRHandler iohA;
  IOHDF5Handler iohB;
  bool is_hdf5 = false;

  string::size_type pos = filename.find("[");
  string basename;
  if (pos != string::npos) {
    basename = filename.substr(0, pos);
  } else {
    basename = filename;
  }
  // Try binary format first
  if (iohA.peekID(ID, basename)) {
    is_hdf5 = false;
  }
  // If that fails, try HDF5 format
  else if (iohB.peekID(ID, basename)) {
    is_hdf5 = true;
  } else {
    m_logger << "Error: could not extract ID string from file " << filename
             << endl;
    throw(std::runtime_error("Could not extract ID string from file"));
  }

  ID = tidyString(ID);
  string sID("Sigmond--SamplingsFile");
  if (ID != sID) {
    m_logger << endl << "This is NOT a Sigmond samplings file" << endl;
    m_logger << "File ID found: <" << ID << ">" << endl;
    m_logger << "Expected ID: <" << sID << ">" << endl;
    throw(std::runtime_error("Invalid Sigmond samplings file"));
  }

  // Now open the file using IOMap which will auto-detect the format
  IOMap<MCObsInfo, Vector<double>> iom;
  iom.openReadOnly(filename, sID);
  XMLHandler xmlh;
  xmlh.set_from_string(iom.getHeader());
  MCBinsInfo bins_info(xmlh);
  MCSamplingInfo samp_info(xmlh, bins_info);
  return make_pair(bins_info, samp_info);
}

// ***************************************************************************************
