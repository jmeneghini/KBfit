#ifndef KB_OBS_HANDLER_H
#define KB_OBS_HANDLER_H
#include "obs_get_handler.h"

// ************************************************************************************
// * *
// *   The class "KBObsHandler", a container class for Monte Carlo non-simple *
// *   observables, is defined in this file.  This is one of the most important
// *
// *   and central classes in "KBfit".  An object of this class is used to get
// the    *
// *   Monte Carlo resamplings, store them in memory, and perform statistical *
// *   analysis on the data.  When calculations are carried out on the data or *
// *   fits are performed, objects of this class store them so that bootstrap
// and     *
// *   jackknife errors can be subsequently computed.  Objects of this class
// also     *
// *   carry out reading and writing of bootstrap/jackknife results to and from
// *
// *   files. *
// * *
// *   An important member of a "KBObsHandler" object include: a ppointer to *
// *   an object of the "MCObsGetHandler" class in the member "m_in_handler" *
// *   which is used for getting the data from files.  A "KBObsHandler" object *
// *   contains the Monte Carlo data in the map: *
// * *
// *     map<KBObsInfo,pair<RVector,uint> > m_samples; *
// * *
// *   Resampling values of non-simple observables are stored in "m_samples". *
// *   In these maps, the keys are of class "KBObsInfo" and the *
// *   values are pairs of RVectors and unsigned integers. The RVectors contain
// *
// *   the resampling values: the 0th element is always the value for the full *
// *   sampling.  For "Jackknife" mode, the index "k" varies from 0 to Nbins,
// the     *
// *   number of bins, where k=0 is the full sample, k=1 is the first jackknife
// *
// *   sampling (bin "k" removed), and so on.  The full sample, stored in the
// k=0     *
// *   element, is used to simplify the computation of each jackknife sample. *
// *   For "Bootstrap" mode, the index "k" ranges from 0 to Nboot, the number of
// *
// *   bootstrap resamplings, where k=0 is the full sample, k=1 is the first *
// *   bootstrap sampling, and so on.  The unsigned integer in each pair
// indicates    *
// *   how many resamplings have already been calculated and stored in memory: *
// *   uncalculated values are stored as quiet NaNs; once this unsigned integer
// *
// *   equals the size of associated RValue, then all samplings are available. *
// * *
// *   The Monte Carlo resamplings are stored in memory, so beware of memory *
// *   exhaustion for a large number of observables.  No routines are provided *
// *   for explicitly reading data, but computing means and covariances causes *
// *   data to be read from file if not already in memory.  A "clearSamplings" *
// *   member is provided to release memory. *
// * *
// *   Data analysis and reading/writing may be accomplished using only a *
// *   **single** "MCSamplingInfo".  Multiple ensembles are allowed, but the *
// *   bootstrap resampling must be the same for all ensembles.  If the sampling
// *
// *   information is of Jackknife mode, then all data must be jackknife *
// *   resamplings, and the jackknife resampling must be the same for all *
// *   observables related to the same ensemble, and the number of bins must be
// *
// *   the same for all ensembles.  The member "setSamplingInfo" can be called
// to     *
// *   change the sampling information, but this clears all memory and
// essentially    *
// *   completely resets the handler. *
// * *
// *   For observables, such as fit parameters, that are not associated with *
// *   any particle Monte Carlo ensemble, a special "MCEnsembleInfo" with *
// *   id string "indep" should be used to indicate an observable that is *
// *   independent of the ensemble. *
// * *
// *   Usage: *
// * *
// *   (1) The constructor requires an "MCSamplingInfo" object.  The constructor
// *
// *       allocates an internal "MCObsGetHandler" object. *
// * *
// *   (2) Information about current parameters can be obtained as below. *
// * *
// *       KH.getNumberOfMeasurements(); *
// *       KH.getNumberOfBins(); *
// *       KH.getNumberOfBootstrapResamplings(); *
// *       KH.getBootstrapper();  // returns const reference *
// * *
// *   (3) Erasing all data corresponding to a particular "KBObsInfo" *
// *   or all data can be accomplished as shown below. *
// * *
// *       KBObsInfo obskey; *
// *       KH.eraseSamplings(obskey);  // erases resamplings for one observable
// *
// *       KH.clearSamplings();        // clears all resamplings *
// * *
// *   (4) Iterating over the resamplings is often needed, either for getting *
// *   the values or putting the values into memory.  Use an integer index. *
// *   Value 0 is for the full sampling. Index 1 is the first resampling, and *
// *   so on. *
// * *
// *   (5) Input/output of all samplings (including full estimates): *
// *   To write samplings (the **default** mode only) that are already in memory
// *
// *   (all samplings include **full* esimates must be available) to file, use *
// * *
// *      set<KBObsInfo> obskeys; .... *
// *      string filename; *
// *      SamplingMode mode=Bootstrap; *
// *      XMLHandler xmlout;   (for output) *
// *      bool overwrite=false;  (file overwriting) *
// *      KH.writeSamplingValuesToFile(obskeys,filename,mode,xmlout,overwrite);
// *
// * *
// *   If "filename" does not exist, it will be created.  Existing files are
// never    *
// *   erased.  New records are added to the files.  If the key of a record to
// be put *
// *   already exists in a file, the put will only occur if "overwrite" is
// specified  *
// *   AND the size of the data to be put does not exceed the size of the data *
// *   already in the file for that key. *
// * *
// *   The format of the files is that for an IOMap.  The header for the *
// *   file contains the information *
// * *
// *       <SigmondSamplingsFile> *
// *          <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo> *
// *          <NumberOfMeasurements>551</NumberOfMeasurements> *
// *          <NumberOfBins>274</NumberOfBins> *
// *          <TweakEnsemble>    (optional) *
// *             <Rebin>2</Rebin> *
// *             <OmitConfig>0</OmitConfig> *
// *             <OmitConfig>1</OmitConfig> *
// *             <OmitConfig>78</OmitConfig> *
// *          </TweakEnsemble> *
// *          <JackKnifeMode/>  or  <Bootstrapper> *
// *                            <NumberResamplings>2048</NumberResamplings> *
// *                                  <Seed>6754</Seed> *
// *                                  <BootSkip>127</BootSkip> *
// *                                </Bootstrapper> *
// *       </SigmondSamplingsFile> *
// * *
// *   The record keys are unsigned integers, and the values contain *
// *   the MCObservable XML string, followed by the vector of samplings. *
// *   An IOMap requires the keys to all have the same size in bytes, so *
// *   an KBObsInfo cannot be used as a key since its encoding can be *
// *   different sizes. *
// * *
// * *
// ************************************************************************************

class KBObsHandler {

  MCObsGetHandler* m_in_handler; // the input handler class
  std::map<KBObsInfo, std::pair<RVector, uint>>
      m_samples; // for storing resamplings
  std::stringstream m_logger;

  // prevent copying
#ifndef NO_CXX11
  KBObsHandler() = delete;
  KBObsHandler(const KBObsHandler& indata) = delete;
  KBObsHandler& operator=(const KBObsHandler& indata) = delete;
#else
  KBObsHandler();
  KBObsHandler(const KBObsHandler& indata);
  KBObsHandler& operator=(const KBObsHandler& indata);
#endif

public:
  KBObsHandler(const MCSamplingInfo& sampinfo);

  ~KBObsHandler();

  void
  setSamplingInfo(const MCSamplingInfo&
                      sampinfo); // if new mode != current, does complete clear

  SamplingMode getSamplingMode() const;

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

  bool isJackknifeMode() const;

  bool isBootstrapMode() const;

  std::set<MCEnsembleInfo> getEnsembleInfos() const;

  const MCBinsInfo& getBinsInfo(const KBObsInfo& obskey) const;

  const MCBinsInfo& getBinsInfo(const MCEnsembleInfo& mcens) const;

  uint getNumberOfBins(const MCEnsembleInfo& mcens) const;

  const MCSamplingInfo& getSamplingInfo() const;

  unsigned int getNumberOfResamplings() const;

  void clearSamplings(); // clears memory, but not underlying MCObsGetHandler;

  void eraseSamplings(const KBObsInfo& obskey);

  void clearGetter(); // clears underlying MCObsGetHandler

  void clear(); // clears memory AND underlying MCObsGetHandler

  bool queryFullAndSamplings(const KBObsInfo& obskey);

  const RVector& getFullAndSamplingValues(const KBObsInfo& obskey);

  // returns null pointer if not available
  const RVector* getFullAndSamplingValuesMaybe(const KBObsInfo& obskey);

  double getFullSampleValue(const KBObsInfo& obskey);

  double getSamplingValue(const KBObsInfo& obskey, int samp_index);

  bool getSamplingValueMaybe(const KBObsInfo& obskey, int samp_index,
                             double& result);

  void putSamplingValue(const KBObsInfo& obskey, int samp_index, double value,
                        bool overwrite = true);

  const RVector& putFullAndSamplings(const KBObsInfo& obskey,
                                     const RVector& samps,
                                     bool overwrite = true);

  void putFixedValue(const KBObsInfo& obskey, const double& fixedvalue,
                     uint nsamplings, bool overwrite = true);
  /*
     RVector getSamplingValues(const KBObsInfo& obskey);

     void getSamplingValues(const KBObsInfo& obskey, RVector& samp);
  */

  MCEstimate getEstimate(const KBObsInfo& obskey);

  double getCovariance(const KBObsInfo& obskey1, const KBObsInfo& obskey2);

  double getStandardError(const KBObsInfo& obskey);

  /*
         // read all samplings from file and put into memory (second version
         // only reads those records matching the KBObsInfo objects in
     "obskeys")

     void readSamplingValuesFromFile(const std::string& filename,
                                     XMLHandler& xmlout);

     void readSamplingValuesFromFile(const std::set<KBObsInfo>& obskeys,
                                     const std::string& filename,
                                     XMLHandler& xmlout);
  */
  // write from memory into file (only if all samplings available)
  // All observables in "obskeys" must be from the same ensemble.

  void writeSamplingValuesToFile(const std::set<KBObsInfo>& obskeys,
                                 const std::string& filename,
                                 XMLHandler& xmlout, bool overwrite = false);

  const std::stringstream& getCurrentLog() const { return m_logger; }

  void clearLog() {
    m_logger.str(std::string());
    m_logger.clear();
  }

private:
  /*
     void assert_simple(const KBObsInfo& obskey, const std::string& name);

     const RVector& get_bins(const KBObsInfo& obskey);
  */
  const std::pair<RVector, uint>*
  get_full_and_sampling_values(const KBObsInfo& obskey);

  /*
     const RVector& put_samplings_in_memory(const KBObsInfo& obskey,
                          const RVector& samplings,
                          std::map<KBObsInfo,std::pair<RVector,uint> >
     *samp_ptr);

     void put_a_sampling_in_memory(const KBObsInfo& obskey, uint sampling_index,
                          double value, bool overwrite, uint sampling_max,
                          std::map<KBObsInfo,std::pair<RVector,uint> >
     *samp_ptr);

     void calc_simple_jack_samples(const RVector& bins, RVector& samplings);

     void calc_simple_boot_samples(const RVector& bins, RVector& samplings);
  */

  void jack_analyze(const RVector& sampvals, MCEstimate& result);

  void boot_analyze(RVector& sampvals, MCEstimate& result);

public:
  double jack_covariance(const RVector& sampvals1, const RVector& sampvals2);

  double boot_covariance(const RVector& sampvals1, const RVector& sampvals2);

  friend class TaskHandler;
};

// *************************************************************
#endif
