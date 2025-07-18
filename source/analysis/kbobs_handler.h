#ifndef KB_OBS_HANDLER_H
#define KB_OBS_HANDLER_H
#include "obs_get_handler.h"

/**
 * @brief Container class for Monte Carlo non-simple observables in KBfit
 *
 * The KBObsHandler class is one of the most important and central classes in
 * KBfit. It is used to get Monte Carlo resamplings, store them in memory, and
 * perform statistical analysis on the data. When calculations are carried out
 * on the data or fits are performed, objects of this class store them so that
 * bootstrap and jackknife errors can be subsequently computed. Objects of this
 * class also carry out reading and writing of bootstrap/jackknife results to
 * and from files.
 *
 * @details
 *
 * **Core Data Storage:**
 *
 * The main data storage is in the member:
 * ```cpp
 * map<KBObsInfo, pair<RVector, uint>> m_samples;
 * ```
 *
 * Resampling values of non-simple observables are stored in `m_samples`, where:
 * - **Keys**: KBObsInfo objects identifying the observable
 * - **Values**: Pairs of RVectors and unsigned integers
 *   - **RVector**: Contains resampling values (0th element is full sampling)
 *   - **uint**: Number of resamplings calculated and stored in memory
 *
 * **Resampling Index Structure:**
 *
 * - **Jackknife mode**: Index k varies from 0 to Nbins
 *   - k=0: Full sample
 *   - k=1: First jackknife sampling (bin 1 removed)
 *   - k=2: Second jackknife sampling (bin 2 removed), etc.
 *
 * - **Bootstrap mode**: Index k ranges from 0 to Nboot
 *   - k=0: Full sample
 *   - k=1: First bootstrap sampling
 *   - k=2: Second bootstrap sampling, etc.
 *
 * **Memory Management:**
 * - Uncalculated values are stored as quiet NaNs
 * - Once the unsigned integer equals the RVector size, all samplings are
 * available
 * - Monte Carlo resamplings are stored in memory - beware of memory exhaustion
 * for large numbers of observables
 * - Use clearSamplings() to release memory when needed
 *
 * **Important Members:**
 * - `m_in_handler`: Pointer to MCObsGetHandler object for getting data from
 * files
 * - `m_samples`: Map storing the resampling data
 *
 * **Sampling Requirements:**
 * - Data analysis uses only a **single** MCSamplingInfo
 * - Multiple ensembles allowed, but bootstrap resampling must be the same for
 * all
 * - For Jackknife mode: all data must be jackknife resamplings with same
 * parameters
 * - Number of bins must be the same for all ensembles
 * - Use setSamplingInfo() to change sampling (clears all memory)
 *
 * **Independent Observables:**
 * For observables like fit parameters not associated with any Monte Carlo
 * ensemble, use a special MCEnsembleInfo with id string "indep".
 *
 * @section usage_examples Usage Examples
 *
 * **1. Constructor Usage:**
 * ```cpp
 * MCSamplingInfo sampling_info;
 * KBObsHandler KH(sampling_info);  // Allocates internal MCObsGetHandler
 * ```
 *
 * **2. Getting Current Parameters:**
 * ```cpp
 * uint num_measurements = KH.getNumberOfMeasurements();
 * uint num_bins = KH.getNumberOfBins();
 * uint num_bootstrap = KH.getNumberOfBootstrapResamplings();
 * const Bootstrapper& bs = KH.getBootstrapper();  // returns const reference
 * ```
 *
 * **3. Memory Management:**
 * ```cpp
 * KBObsInfo obskey;
 * KH.eraseSamplings(obskey);  // erases resamplings for one observable
 * KH.clearSamplings();        // clears all resamplings
 * ```
 *
 * **4. Resampling Iteration:**
 * Use integer index where:
 * - Value 0: Full sampling
 * - Index 1: First resampling
 * - Index 2: Second resampling, etc.
 *
 * **5. File I/O Operations:**
 * ```cpp
 * set<KBObsInfo> obskeys;
 * string filename;
 * SamplingMode mode = Bootstrap;
 * XMLHandler xmlout;
 * bool overwrite = false;
 * KH.writeSamplingValuesToFile(obskeys, filename, mode, xmlout, overwrite);
 * ```
 *
 * **File Behavior:**
 * - Non-existent files are created
 * - Existing files are never erased
 * - New records are added to files
 * - Existing records are only overwritten if `overwrite=true` AND new data size
 * ≤ existing data size
 *
 * @section file_format File Format
 *
 * Files use IOMap format with header containing:
 *
 * @code{.xml}
 * <SigmondSamplingsFile>
 *   <MCEnsembleInfo>clover_s24_t128_ud840_s743</MCEnsembleInfo>
 *   <NumberOfMeasurements>551</NumberOfMeasurements>
 *   <NumberOfBins>274</NumberOfBins>
 *   <TweakEnsemble>    <!-- optional -->
 *     <Rebin>2</Rebin>
 *     <OmitConfig>0</OmitConfig>
 *     <OmitConfig>1</OmitConfig>
 *     <OmitConfig>78</OmitConfig>
 *   </TweakEnsemble>
 *   <JackKnifeMode/>  <!-- or -->
 *   <Bootstrapper>
 *     <NumberResamplings>2048</NumberResamplings>
 *     <Seed>6754</Seed>
 *     <BootSkip>127</BootSkip>
 *   </Bootstrapper>
 * </SigmondSamplingsFile>
 * @endcode
 *
 * **Record Structure:**
 * - **Keys**: Unsigned integers (IOMap requires same-size keys)
 * - **Values**: MCObservable XML string + vector of samplings
 *
 * @note KBObsInfo cannot be used directly as keys since its encoding can vary
 * in size
 *
 * @warning Computing means and covariances causes data to be read from file if
 * not already in memory
 *
 * @warning Be careful of memory exhaustion when handling large numbers of
 * observables
 */

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
  double get_average(const RVector& sampvals);

  double jack_covariance(const RVector& sampvals1, const RVector& sampvals2);

  double boot_covariance(const RVector& sampvals1, const RVector& sampvals2);

  friend class TaskHandler;
};

// *************************************************************
#endif
