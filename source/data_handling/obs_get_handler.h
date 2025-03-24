#ifndef OBS_GET_HANDLER_H
#define OBS_GET_HANDLER_H

#include "ensemble_info.h"
#include "kbobs_info.h"
#include "sampling_info.h"
#include "samplings_handler.h"

// *******************************************************************
// *                                                                 *
// *  "MCObsGetHandler" handles access to all Monte Carlo            *
// *  observables (fit parameters) stored in files.  It is a         *
// *  mid-level routine: it is meant to be used by the higher level  *
// *  "KBObsHandler" and it uses the lower-level                     *
// *  "SamplingsGetHandler" handler.  This handler keeps track       *
// *  of which Monte Carlo observables are available and reads       *
// *  them from file when requested.  Results are not stored in      *
// *  memory: the higher level "KBObsHandler" is responsible for     *
// *  storing samplings in memory.                                   *
// *                                                                 *
// *  This handler is very different from its "SigMonD" counterpart. *
// *  This handler only deals with samplings files, but both         *
// *  bootstrap and jackknife resampling are allowed.  When          *
// *  constructed, on object of class "MCObsGetHandler" must be      *
// *  set to either "bootstrap" or "jackknife" resampling mode.      *
// *  All files must match the sampling info.  Since we wish to use  *
// *  multiple ensembles, all bootstrap resampling parameters        *
// *  must be the same for all observables in order to perform the   *
// *  statistical analysis.  For jackknife resampling, the parameters*
// *  must all be the same for EACH ensemble.  Covariances between   *
// *  observables on different ensembles are taken to be zero.       *
// *                                                                 *
// *  The most important member routines are                         *
// *                                                                 *
// *     void getSamplings(const MCObsInfo& obsinfo, RVector& samp); *
// *     bool querySamplings(const MCObsInfo& obsinfo);              *
// *                                                                 *
// *  The "get" routines above throw an exception if the observable  *
// *  is not available.                                              *
// *                                                                 *
// *  The constructor takes an XMLHandler as input parameter, as     *
// *  well as an "MCBinsInfo" object and an "MCSamplingInfo"         *
// *  object.  The rebinning and omission of configurations in       *
// *  "MCBinsInfo" must be consistent the header information in the  *
// *  files: same omissions, but rebin factor can be integer         *
// *  multiple of rebin factor in file.                              *
// *                                                                 *
// *  Similarly for sampling data, the precise resampling used       *
// *  as specified in the "MCSamplingInfo" object must match the     *
// *  header information in the sampling files.  Only one sampling   *
// *  mode (Jackknife or Bootstrap) can be handled.  The XML needed  *
// *  for the constructor must have the form                         *
// *                                                                 *
// *     <KBObservables>                                             *
// *                                                                 *
// *       <MCSamplingInfo>...</MCSamplingInfo>                      *
// *                                                                 *
// *       <SamplingData>                                            *
// *          <FileName>...</FileName>                               *
// *             ....                                                *
// *       </SamplingData>                                           *
// *                                                                 *
// *                                                                 *
// *       <Specifications>                                          *
// *             ...                                                 *
// *           specifications of observables (optional)              *
// *             ...                                                 *
// *       </Specifications>                                         *
// *                                                                 *
// *       <Verbose/>  (optional)                                    *
// *                                                                 *
// *     </KBObservables>                                            *
// *                                                                 *
// *  If no observables are specified, then all observables found    *
// *  while reading the files will be used.  If any observables      *
// *  ARE specified, then only those observables will be             *
// *  considered for input: files corresponding to other             *
// *  observables will be ignored, and an exception is thrown if     *
// *  the file containing the data for any requested observable      *
// *  cannot be found.                                               *
// *                                                                 *
// *  Observables can be specified inside a <Specifications>         *
// *  tag as follows:                                                *
// *                                                                 *
// *       <ObsSamplings>                                            *
// *           <MCEnsembleInfo>...</MCEnsembleInfo>                  *
// *           <MCObservable>...</MCObservable>                      *
// *           <MCObservable>...</MCObservable>                      *
// *                ...                                              *
// *       </ObsSamplings>                                           *
// *       <ObsSamplings> ...</ObsSamplings>  (for each ensemble)    *
// *       <ObsSamplings>                                            *
// *           <MCEnsembleInfo>...</MCEnsembleInfo>                  *
// *           <All/>  (all keys for this ensemble)                  *
// *       </ObsSamplings>                                           *
// *                                                                 *
// *                                                                 *
// *  Sampling files can be added and removed later using the        *
// *  "connect" and "disconnect" members below.                      *
// *                                                                 *
// *                                                                 *
// ******************************************************************

// *******************************************************************
// *                                                                 *
// *  Implementation notes:                                          *
// *                                                                 *
// *  An object of this "MCObsGetHandler" class maintains a map with *
// *  record keys of class "MCEnsembleInfo" and values of type       *
// *  pointer to a "SamplingsGetHandler".  The objects pointed to    *
// *  perform the actual data reading from files but they do NOT     *
// *  store data in memory.  "MCObsGetHandler" objects also do       *
// *  NOT store any data in memory.  Instead, the higher level       *
// *  "KBObsHandler" object maintains the input data in memory, in a *
// *  format that facilitates evaluating mean values, covariances,   *
// *  bootstrap errors, and so on.  The class "MCObsGetHandler" is   *
// *  meant to be used as a data handler for an object of the        *
// *  "KBObsHandler" class, with "MCEnsembleInfo" and "MCObsInfo"    *
// *  as the record key type.                                        *
// *                                                                 *
// *                                                                 *
// *******************************************************************

class MCObsGetHandler {

  std::map<MCEnsembleInfo, SamplingsGetHandler*> m_sampsdh;
  MCSamplingInfo* m_sampinfo;
  std::map<MCEnsembleInfo, MCBinsInfo> m_bininfos;
  SamplingMode m_sampmode;
  std::stringstream m_logger;

  // Prevent copying ... handler might contain large
  // amounts of data

#ifndef NO_CXX11
  MCObsGetHandler() = delete;
  MCObsGetHandler(const MCObsGetHandler&) = delete;
  MCObsGetHandler& operator=(const MCObsGetHandler&) = delete;
#else
  MCObsGetHandler();
  MCObsGetHandler(const MCObsGetHandler&);
  MCObsGetHandler& operator=(const MCObsGetHandler&);
#endif

public:
  MCObsGetHandler(XMLHandler& xml_in);

  MCObsGetHandler(const MCSamplingInfo& sampinfo);

  MCObsGetHandler(const MCSamplingInfo& sampinfo,
                  const std::set<std::string>& sampfiles, bool verbose = false);

  MCObsGetHandler(
      const MCSamplingInfo& sampinfo, const std::set<std::string>& sampfiles,
      const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep,
      bool verbose = false);

  ~MCObsGetHandler();

  void connectSamplingFiles(const std::set<std::string>& sampfiles,
                            bool verbose = false);

  // if map "keys_to_keep" is empty:
  //    -- keep all keys encountered for ALL ensembles encountered
  // else if "keys_to_keep" DOES have an entry for an ensemble:
  //    -- if entry is empty set, keep all keys for that ensemble
  //    -- if entry is nonempty set, keep only those keys in the set
  // else if "keys_to_keep" does NOT have an entry for an ensemble:
  //    -- no keys are to be kept for that ensemble

  void connectSamplingFiles(
      const std::set<std::string>& sampfiles,
      const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep,
      bool verbose = false);

  void disconnectSamplingFiles(const std::set<std::string>& sampfiles);

  void disconnectAllSamplingFiles();

  // "setSamplingInfo" below does a complete clear if the
  // new mode differs from current; nothing happens if new=current mode

  void setSamplingInfo(const MCSamplingInfo& sampinfo);

  const MCBinsInfo& getBinsInfo(const MCEnsembleInfo& mcens) const;

  const MCSamplingInfo& getSamplingInfo() const;

  uint getNumberOfResamplings() const;

  void close();

  void getFileMap(XMLHandler& xmlout) const;

  SamplingMode getSamplingMode() const;

  bool isJackknifeMode() const;

  bool isBootstrapMode() const;

  std::set<MCEnsembleInfo> getEnsembleInfos() const;

  void getSamplings(const KBObsInfo& obsinfo, RVector& samp) const;

  bool getSamplingsMaybe(const KBObsInfo& obsinfo, RVector& samp) const;

  bool querySamplings(const KBObsInfo& obsinfo) const;

  void clear();

  const std::stringstream& getCurrentLog() const { return m_logger; }

  void clearLog() {
    m_logger.str(std::string());
    m_logger.clear();
  }

private:
  std::pair<MCBinsInfo, MCSamplingInfo>
  get_info_from_file(const std::string& filename);

  void initialize(const MCSamplingInfo& sampinfo);

  void connect_a_samplings_file(
      const std::string& file_name,
      const std::map<MCEnsembleInfo, std::set<MCObsInfo>>& keys_to_keep);
};

// ***************************************************************
#endif
