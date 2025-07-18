#include "chisq_detres.h"
#include "chisq_fit.h"
#include "root_finder.h"
#include "task_handler.h"
#include <filesystem>
#include <iomanip>

using namespace std;

/**
 * @brief Prints the Omega function for Elab between minimum and maximum values at fixed increments
 * 
 * This function evaluates and outputs the Omega function (or determinant) for lab-frame energies
 * between user-specified minimum and maximum values. Output for each KBBlock is sent to a separate 
 * file with a common stub and different integer suffixes.
 * 
 * @param xmltask XML task handler containing input parameters
 * @param xmlout XML output handler for results
 * @param taskcount Task counter for file naming
 * 
 * @details
 * 
 * **Output Formats:**
 * - `OutputMode = "full"`: Format is `Elab/mref  value_from_full_sample`
 * - `OutputMode = "resampled"`: Format is `Elab/mref  value_from_full_sample upper_error down_error`
 * 
 * **Requirements:**
 * - An L³ spatial lattice is required
 * - Reference scale times temporal lattice spacing must be specified in `<ReferenceMassTimeSpacingProduct>`
 * - Either `<KtildeMatrixInverse>` or `<KtildeMatrix>` must be provided
 * 
 * **Configuration Options:**
 * - If `<OmegaMu>` is specified, the Omega function is used; otherwise, the determinant itself is used
 * - For anisotropic lattices, specify `<LatticeAnisotropy>` (spatial/temporal spacing ratio)
 * - Lab-frame energies and masses can be input as reference mass ratios or time spacing products
 * 
 * @note Energy format should be specified in `<DefaultEnergyFormat>` as either 
 *       "reference_ratio" or "time_spacing_product"
 * 
 * **XML Input Format:**
 * 
 * @code{.xml}
 * <Task>
 *   <Action>DoPrint</Action>
 *   <OutputMode>full</OutputMode> <!-- or "resampled" -->
 *   <OmegaMu>8.0</OmegaMu> <!-- optional -->
 *   <QuantizationCondition>StildeCB</QuantizationCondition>
 *   <!-- Options: StildeCB, StildeinvCB, KtildeB, KtildeinvB -->
 *   
 *   <RootFinder>
 *     <LabFrameEnergyMin>1.10</LabFrameEnergyMin>
 *     <LabFrameEnergyMax>2.30</LabFrameEnergyMax>
 *     
 *     <AdaptiveBracket>
 *       <!-- All tags below are OPTIONAL -->
 *       <InitialStepSize>1e-2</InitialStepSize>
 *       <AbsXTolerance>1e-6</AbsXTolerance>
 *       <AbsResidualTolerance>1e-10</AbsResidualTolerance>
 *       <MinStepSize>1e-4</MinStepSize>
 *       <MaxStepSize>5e-2</MaxStepSize>
 *       <StepScaleLimit>3.0</StepScaleLimit>
 *       <PlateauMod2Threshold>1e-4</PlateauMod2Threshold>
 *       <PlateauCountBeforeJump>4</PlateauCountBeforeJump>
 *     </AdaptiveBracket>
 *   </RootFinder>
 *   
 *   <KtildeMatrixInverse> <!-- or <KtildeMatrix> -->
 *     <!-- K-matrix specification... -->
 *   </KtildeMatrixInverse>
 *   
 *   <DefaultEnergyFormat>reference_ratio</DefaultEnergyFormat>
 *   <!-- or time_spacing_product -->
 *   
 *   <MCEnsembleParameters>...</MCEnsembleParameters>
 *   <!-- One for each Monte Carlo ensemble -->
 *   
 *   <KBBlock>...</KBBlock>
 *   <!-- One for each KB quantization block -->
 *   
 *   <KBObservables>
 *     <MCSamplingInfo>...</MCSamplingInfo>
 *     <SamplingData>
 *       <FileName>...</FileName>
 *     </SamplingData>
 *   </KBObservables>
 * </Task>
 * @endcode
 * 
 * **AdaptiveBracket Parameter Mapping:**
 * | XML Tag                      | AdaptiveBracketConfig Field |
 * |------------------------------|----------------------------|
 * | `<InitialStepSize>`          | `initial_step_size`        |
 * | `<AbsXTolerance>`            | `x_tol`                    |
 * | `<AbsResidualTolerance>`     | `zero_tol`                 |
 * | `<MinStepSize>`              | `min_step_size`            |
 * | `<MaxStepSize>`              | `max_step_size`            |
 * | `<StepScaleLimit>`           | `step_scale_limit`         |
 * | `<PlateauMod2Threshold>`     | `plateau_mod2_threshold`   |
 * | `<PlateauCountBeforeJump>`   | `plateau_count_before_jump`|
 * 
 * **Monte Carlo Ensemble Configuration:**
 * @code{.xml}
 * <MCEnsembleParameters>
 *   <MCEnsembleInfo>...</MCEnsembleInfo>
 *   <ReferenceMassTimeSpacingProduct>
 *     <MCObs>...</MCObs>
 *   </ReferenceMassTimeSpacingProduct>
 *   <LatticeAnisotropy> <!-- optional: a_s/a_t, unity if absent -->
 *     <MCObs>...</MCObs>
 *   </LatticeAnisotropy>
 *   <ParticleMass>
 *     <Name>pion</Name> <!-- should match names used in K-matrix -->
 *     <MCObs>...</MCObs>
 *   </ParticleMass>
 *   <!-- Additional particle masses... -->
 * </MCEnsembleParameters>
 * @endcode
 * 
 * **MCObsInfo Specification:**
 * Can use either short or long form (must be nonsimple and real):
 * @code{.xml}
 * <!-- Long form -->
 * <MCObservable>
 *   <ObsName>T1up_Energy</ObsName> <!-- 32 chars max, no blanks -->
 *   <Index>3</Index> <!-- optional nonneg integer, default 0 -->
 * </MCObservable>
 * 
 * <!-- Short form -->
 * <MCObs>T1up_Energy 3</MCObs>
 * @endcode
 * 
 * **KB Quantization Block Configuration:**
 * @code{.xml}
 * <KBBlock>
 *   <MCEnsembleInfo>...</MCEnsembleInfo>
 *   <BoxQuantization>
 *     <TotalMomentumRay>ar</TotalMomentumRay>
 *     <TotalMomentumIntSquared>0</TotalMomentumIntSquared>
 *     <LGIrrep>T1u</LGIrrep>
 *     <LmaxValues>5 3</LmaxValues> <!-- one for each decay channel -->
 *   </BoxQuantization>
 *   <LabFrameEnergyMin>...</LabFrameEnergyMin>
 *   <LabFrameEnergyMax>...</LabFrameEnergyMax>
 *   <LabFrameEnergyInc>...</LabFrameEnergyInc>
 * </KBBlock>
 * @endcode
 * 
 * **K-Matrix Specification:**
 * @code{.xml}
 * <KtildeMatrixInverse> <!-- or <KtildeMatrix> -->
 *   <DecayChannels>
 *     <DecayChannelInfo>
 *       <Particle1Name>pion</Particle1Name>
 *       <Spin1TimesTwo>0</Spin1TimesTwo>
 *       <Identical/> <!-- if identical, omit tags below -->
 *       <Particle2Name>eta</Particle2Name>
 *       <Spin2TimesTwo>2</Spin2TimesTwo>
 *       <IntrinsicParities>same</IntrinsicParities> <!-- or "opposite" -->
 *     </DecayChannelInfo>
 *     <!-- Additional decay channels... -->
 *   </DecayChannels>
 *   
 *   <Element>
 *     <KElementInfo>
 *       <JTimesTwo>2</JTimesTwo>
 *       <KIndex>L(1) 2S(0) chan(0)</KIndex>
 *       <KIndex>L(1) 2S(0) chan(0)</KIndex>
 *     </KElementInfo>
 *     <FitForm>
 *       <Polynomial><Powers>1 3</Powers></Polynomial>
 *     </FitForm>
 *   </Element>
 *   <!-- Additional elements... -->
 *   
 *   <StartingValues>
 *     <KFitParamInfo>
 *       <PolynomialTerm>
 *         <Power>3</Power>
 *         <KElementInfo>...</KElementInfo>
 *       </PolynomialTerm>
 *       <StartingValue>0.4534</StartingValue>
 *     </KFitParamInfo>
 *     <!-- Additional starting values... -->
 *   </StartingValues>
 * </KtildeMatrixInverse>
 * @endcode
 * 
 * @warning The reference mass, particle masses, and anisotropy for each ensemble can be set 
 *          to fixed values by replacing `<MCObs>`/`<MCObservable>` tags with 
 *          `<FixedValue>1.1123</FixedValue>` tags. However, this practice can lead to 
 *          erroneous error estimates and should be used with caution.
 * 
 * @note Order matters in decay channel specification: first `<DecayChannelInfo>` is channel 0,
 *       second is channel 1, etc. In the K-matrix, channels are referenced by index 0, 1, ...
 * 
 * @note Fixed values for particle masses always refer to energy ratios over the reference mass.
 */

void TaskHandler::doPrint(XMLHandler& xmltask, XMLHandler& xmlout,
                          int taskcount) {
  try {
    xmlout.set_root("PrintOmega");
    stringstream logger;
    logger.precision(12);
    KtildeMatrixCalculator* Kmat = 0;
    KtildeInverseCalculator* Kinv = 0;

    // check if we're printing eigenvalues
    int count_print_eigenvals =
        xmltask.count_among_children("PrintEigenvalues");
    if (count_print_eigenvals > 1) {
      throw(std::invalid_argument(
          "Multiple PrintEigenvalues tags cannot be present"));
    }
    bool do_print_eigenvals = count_print_eigenvals;

    string outmode = "full";
    xmlreadif(xmltask, "OutputMode", outmode, "doPrint");
    if ((outmode != "full") && (outmode != "resampled"))
      throw(std::runtime_error("Invalid output mode specified"));

    //  get sampling mode and names of input files
    XMLHandler xmlr(xmltask, "KBObservables");
    MCSamplingInfo sampinfo(xmlr);
    if (sampinfo != m_obs->getSamplingInfo()) {
      throw(std::invalid_argument("KBObservables MCSamplingInfo does not match "
                                  "that of the KBObsHandler"));
    }

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
      logger << "No sampling files given" << endl;
      return;
    }

    //  first set up the K or K^(-1) matrix
    // uint numfitparams;
    int k1 = xmltask.count_among_children("KtildeMatrix");
    int k2 = xmltask.count_among_children("KtildeMatrixInverse");
    if ((k1 + k2) != 1)
      throw(std::invalid_argument(
          "A single KtildeMatrix or KtildeMatrixInverse tag must be present"));
    const vector<DecayChannelInfo>* dcptr;
    if (k1 == 1) {
      Kmat = new KtildeMatrixCalculator(xmltask, true);
      // numfitparams=Kmat->getNumberOfParameters();
      dcptr = &(Kmat->getDecayChannelInfos());
    } else {
      Kinv = new KtildeInverseCalculator(xmltask, true);
      // numfitparams=Kinv->getNumberOfParameters();
      dcptr = &(Kinv->getDecayChannelInfos());
    }

    // get set of particle names
    set<string> pnames;
    for (vector<DecayChannelInfo>::const_iterator dct = dcptr->begin();
         dct != dcptr->end(); dct++) {
      pnames.insert(dct->getParticle1Name());
      pnames.insert(dct->getParticle2Name());
    }

    //  assign mu  (negative value means use determinant itself)
    double omega_mu = -1.0;
    xmlreadif(xmltask, "OmegaMu", omega_mu, "DoPrint");
    logger << "Omega mu = " << omega_mu << endl;

    //  get format for energies/mass (in terms of product with
    //  lattice time spacing or a ratio with reference mass)
    bool energy_ratios = true;
    {
      string reply;
      xmlreadif(xmltask, "DefaultEnergyFormat", reply, "DoPrint");
      if ((reply == "reference_ratio") || (reply.empty()))
        energy_ratios = true;
      else if (reply == "time_spacing_product")
        energy_ratios = false;
      else
        throw(std::invalid_argument("Invalid <DefaultEnergyFormat> tag"));
    }

    // Get root finding info
    int root_counts = xmltask.count_among_children("RootFinder");
    if (root_counts > 1)
      throw(
          std::invalid_argument("Multiple RootFinder tags cannot be present"));
    bool do_root_find = root_counts;
    AdaptiveBracketConfig root_config;
    if (do_root_find) {
      XMLHandler xmlroot(xmltask, "RootFinder");
      root_config = AdaptiveBracketRootFinder::makeConfigFromXML(xmlroot);
    }
    logger << "RootFinder configuration: " << endl;
    logger << root_config.toString() << endl;

    //  Loop over the different MC ensembles to get information
    //  about all needed parameters.  Particle masses, anisotropy,
    //  parameters should be associated with the ensemble.

    map<MCEnsembleInfo, set<MCObsInfo>> needed_keys;
    set<MCEnsembleInfo> ensembles;
    map<MCEnsembleInfo, MCObsInfo> ref_at_mass;
    map<MCEnsembleInfo, MCObsInfo> anisotropy;
    map<MCEnsembleInfo, map<string, MCObsInfo>> particle_masses;
    map<KBObsInfo, double> fixed_values;

    list<XMLHandler> xmlen = xmltask.find("MCEnsembleParameters");
    for (list<XMLHandler>::iterator it = xmlen.begin(); it != xmlen.end();
         ++it) {
      MCEnsembleInfo mcens(*it);
      if (ensembles.find(mcens) != ensembles.end())
        throw(std::invalid_argument("Duplicate MC ensemble"));
      ensembles.insert(mcens);
      set<MCObsInfo>& kset = needed_keys[mcens];
      string pname;
      MCObsInfo rkey;
      // get lattice spacing times reference scale info
      ChiSquare::read_obs(*it, "ReferenceMassTimeSpacingProduct", false, rkey,
                          kset, pname, mcens, fixed_values);
      ref_at_mass.insert(make_pair(mcens, rkey));
      if (!energy_ratios) {
        kset.erase(rkey);
      } // remove ref_at_mass from kset
      // get anisotropy info
      if (xml_tag_count(*it, "LatticeAnisotropy") == 1) {
        ChiSquare::read_obs(*it, "LatticeAnisotropy", false, rkey, kset, pname,
                            mcens, fixed_values);
        anisotropy.insert(make_pair(mcens, rkey));
      }
      // get particle mass infos
      map<string, MCObsInfo> pmap;
      list<XMLHandler> xmlp = it->find("ParticleMass");
      for (list<XMLHandler>::iterator pt = xmlp.begin(); pt != xmlp.end();
           ++pt) {
        ChiSquare::read_obs(*pt, "ParticleMass", true, rkey, kset, pname, mcens,
                            fixed_values);
        if (pmap.find(pname) != pmap.end())
          throw(std::invalid_argument("Duplicate particle masses"));
        pmap.insert(make_pair(pname, rkey));
      }
      //  now check that all decay particles have available masses
      bool pcheck = (pmap.size() >= pnames.size());
      set<string>::const_iterator pnt = pnames.begin();
      while ((pnt != pnames.end()) && (pcheck)) {
        pcheck = (pmap.find(*pnt) != pmap.end());
        ++pnt;
      }
      if (!pcheck)
        throw(std::runtime_error("A particle does not have available mass"));
      particle_masses.insert(make_pair(mcens, pmap));
    }
    if (ensembles.empty())
      throw(std::invalid_argument("No ensembles listed in input XML"));

    //  Loop over the KB quantization blocks to get the lab-frame
    //  energies infos.

    vector<double> elabs_min, elabs_max, elabs_inc;
    vector<MCEnsembleInfo> blockens;
    map<MCEnsembleInfo, uint> ensemble_idmap;
    vector<BoxQuantization*> BQ;

    list<XMLHandler> xmlkb = xmltask.find("KBBlock");
    uint blockcount = 0;
    uint ensemblecount = 0;
    for (list<XMLHandler>::iterator it = xmlkb.begin(); it != xmlkb.end();
         ++it) {
      BoxQuantization* bqptr = new BoxQuantization(*it, Kmat, Kinv);
      BQ.push_back(bqptr);
      MCEnsembleInfo mcens(*it);
      if (ensembles.find(*it) == ensembles.end())
        throw(
            std::invalid_argument("KBBlock associated with unknown ensemble"));
      blockens.push_back(mcens);
      if (ensemble_idmap.find(mcens) == ensemble_idmap.end()) {
        ensemble_idmap.insert(make_pair(mcens, ensemblecount));
        ensemblecount++;
      }
      double elabmin, elabmax, elabinc;
      xmlread(*it, "LabFrameEnergyMin", elabmin, "DoPrint");
      xmlread(*it, "LabFrameEnergyMax", elabmax, "DoPrint");
      xmlread(*it, "LabFrameEnergyInc", elabinc, "DoPrint");
      if ((elabmin >= elabmax) || (elabinc <= 0.0)) {
        throw(std::runtime_error(
            "Invalid energy range specification in DoPrint"));
      }
      elabs_min.push_back(elabmin);
      elabs_max.push_back(elabmax);
      elabs_inc.push_back(elabinc);
      blockcount++;
    }
    if (blockcount == 0) {
      throw(std::runtime_error("No data to analyze"));
    }

    // get QuantizationCondition
    string qctype;
    xmlreadif(xmltask, "QuantizationCondition", qctype, "DoPrint");
    if (qctype.empty()) {
      throw(std::invalid_argument("QuantizationCondition tag must be present"));
    }
    BoxQuantization* bqptr_dummy = BQ[0];
    BoxQuantization::QuantCondType qctype_enum;
    try {
      qctype_enum = bqptr_dummy->getQuantCondTypeFromString(qctype).value();
      if (qctype_enum == BoxQuantization::KtildeB) {
        if (k2 == 1) {
          throw(std::invalid_argument(
              "KtildeMatrixInverse cannot be used with StildeCB or KtildeB"));
        }
      }
      if (qctype_enum == BoxQuantization::KtildeinvB) {
        if (k1 == 1) {
          throw(std::invalid_argument(
              "KtildeMatrix cannot be used with StildeinvCB or KtildeinvB"));
        }
      }
    } catch (std::bad_optional_access&) {
      throw(std::invalid_argument("Invalid QuantizationCondition tag"));
    }

    //  connect files for input

    m_obs->connectSamplingFiles(sampfiles, needed_keys, true);
    logger << m_obs->getCurrentLog().str();
    m_obs->clearLog();

    //  get number of resamplings.  For bootstrap, all ensembles
    //  must have same bootstrap parameters.  For jackknife,
    //  all ensembles must have the SAME number of bins.

    uint nsamplings = m_obs->getNumberOfResamplings();
    if (m_obs->isJackknifeMode()) {
      for (set<MCEnsembleInfo>::const_iterator et = ensembles.begin();
           et != ensembles.end(); ++et) {
        if (m_obs->getNumberOfBins(*et) != nsamplings) {
          throw(std::invalid_argument(
              "All ensembles must have same number of jackknife bins"));
        }
      }
    }

    if (nsamplings == 0)
      throw(std::runtime_error(
          "Samplings numbers do not match for all ensembles"));

    //  insert all of the fixed values into m_obs

    for (map<KBObsInfo, double>::const_iterator fx = fixed_values.begin();
         fx != fixed_values.end(); ++fx) {
      m_obs->putFixedValue(fx->first, fx->second,
                           m_obs->getNumberOfResamplings(), true);
    }

    //  read the reference mass time-spacing products for each ensemble into
    //  memory; evaluate reference lengths and particle masses for each ensemble

    for (map<MCEnsembleInfo, MCObsInfo>::iterator et = ref_at_mass.begin();
         et != ref_at_mass.end(); ++et) {
      const MCEnsembleInfo& mcens = et->first;
      uint nsamp = m_obs->getNumberOfResamplings();
      KBObsInfo atrefmasskey(mcens, et->second);
      const RVector& atrefmass0 = m_obs->getFullAndSamplingValues(atrefmasskey);
      if (atrefmass0.size() != (nsamp + 1))
        throw(std::runtime_error("Resampling size mismatch in KBfit"));
      KBObsInfo scalekey(mcens, MCObsInfo("KBScale"));
      m_obs->putFullAndSamplings(scalekey, atrefmass0, true);
      m_obs->eraseSamplings(atrefmasskey);

      // get the reference length

      uint Llat = mcens.getLatticeXExtent();
      if ((mcens.getLatticeYExtent() != Llat) ||
          (mcens.getLatticeZExtent() != Llat))
        throw(
            std::runtime_error("KBfit only works for LxLxL spatial lattices"));
      const RVector& atrefmass = m_obs->getFullAndSamplingValues(scalekey);
      if (atrefmass.size() != (nsamp + 1))
        throw(std::runtime_error("Resampling size mismatch in KBfit"));
      MCObsInfo obskey("LengthReference");
      KBObsInfo lengthkey(mcens, obskey);
      {
        RVector buff(nsamp + 1);
        if (anisotropy.find(mcens) == anisotropy.end()) { // isotropic case
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat);
        } else {
          KBObsInfo anisotropykey(mcens, anisotropy[mcens]);
          const RVector& anisotropy =
              m_obs->getFullAndSamplingValues(anisotropykey);
          if (anisotropy.size() != (nsamp + 1))
            throw(std::runtime_error("Resampling size mismatch in KBfit"));
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat) * anisotropy[k];
        }
        m_obs->putFullAndSamplings(lengthkey, buff);
      }

      //  get particle masses (form ratios if energy_ratio false)

      map<string, MCObsInfo>& pmap = particle_masses[mcens];
      for (set<string>::const_iterator pt = pnames.begin(); pt != pnames.end();
           ++pt) {
        KBObsInfo masskey(mcens, pmap[*pt]);
        const RVector& mass = m_obs->getFullAndSamplingValues(
            masskey); // reads from file, gets into memory
        bool not_fixed = (fixed_values.find(masskey) == fixed_values.end());
        if ((!energy_ratios) && (not_fixed)) {
          RVector massratio(mass);
          if (massratio.size() != atrefmass.size())
            throw(
                std::runtime_error("Size mismatch while forming mass ratios"));
          for (uint kk = 0; kk < massratio.size(); ++kk)
            massratio[kk] /= atrefmass[kk];
          m_obs->putFullAndSamplings(masskey, massratio, true);
        }
      }
    }

    //  Create output directory structure using the new utility function
    filesystem::path output_path =
        createKBOutputDirectory(m_output_directory, qctype);

    for (uint blocknum = 0; blocknum < BQ.size(); ++blocknum) {
      const MCEnsembleInfo& mcens = blockens[blocknum];
      string omega_filename = "Omega_Block" + make_string(blocknum) + ".csv";
      string eigenvals_filename =
          "Eigenvals_Block" + make_string(blocknum) + ".csv";
      string nis_filename =
          "NonInteractingEnergies_Block" + make_string(blocknum) + ".csv";
      string roots_filename = "Roots_Block" + make_string(blocknum) + ".csv";

      // Create full file paths
      filesystem::path omega_filepath = output_path / omega_filename;
      filesystem::path eigenvals_filepath = output_path / eigenvals_filename;
      filesystem::path nis_filepath = output_path / nis_filename;
      filesystem::path roots_filepath = output_path / roots_filename;

      BoxQuantization* bqptr = BQ[blocknum];
      logger << "Output directory = " << output_path << endl;
      logger << "Filename = " << omega_filename << endl;
      if (do_print_eigenvals) {
        logger << "Filename = " << eigenvals_filename << endl;
      }
      if (do_root_find) {
        logger << "Filename = " << roots_filename << endl;
      }
      logger << "Filename = " << nis_filename << endl;
      logger << mcens.str() << endl;
      logger << "MomRay " << bqptr->getMomRay()
             << "   P^2 = " << bqptr->getTotalMomentumIntegerSquared()
             << " Box Irrep " << bqptr->getLittleGroupBoxIrrep() << endl;
      uint nsamp = m_obs->getNumberOfResamplings();
      map<string, MCObsInfo>& pmap = particle_masses[mcens];

      KBObsInfo lengthkey(mcens, MCObsInfo("LengthReference"));
      const RVector& mrefL = m_obs->getFullAndSamplingValues(lengthkey);
      KBObsInfo scalekey(mcens, MCObsInfo("KBScale"));
      // const RVector& atrefmass=m_obs->getFullAndSamplingValues(scalekey);

      uint nchan = bqptr->getNumberOfDecayChannels();
      vector<const RVector*> particlemass1(nchan), particlemass2(nchan);
      for (uint ci = 0; ci < nchan; ++ci) {
        const DecayChannelInfo& chan = bqptr->getDecayChannelInfo(ci);
        const string& pname1 = chan.getParticle1Name();
        const string& pname2 = chan.getParticle2Name();
        KBObsInfo mass1key(mcens, pmap[pname1]);
        KBObsInfo mass2key(mcens, pmap[pname2]);
        particlemass1[ci] = &(m_obs->getFullAndSamplingValues(mass1key));
        particlemass2[ci] = &(m_obs->getFullAndSamplingValues(mass2key));
      }

      for (int b = 0; b <= nsamp; ++b) { // set masses
        bqptr->setRefMassL(mrefL[b]);
        for (uint ci = 0; ci < bqptr->getNumberOfDecayChannels(); ++ci) {
          bqptr->setMassesOverRef(ci, (*(particlemass1[ci]))[b],
                                  (*(particlemass2[ci]))[b]);
        }
      }

      double emin = elabs_min[blocknum];
      double emax = elabs_max[blocknum];

      // Non interacting energies first for resolution of elabsvals

      string header = "#" + mcens.str() + " # MomRay " + bqptr->getMomRay() +
                      " # P^2 = " +
                      std::to_string(bqptr->getTotalMomentumIntegerSquared()) +
                      " # Box Irrep " + bqptr->getLittleGroupBoxIrrep();

      ofstream fout_nis(nis_filepath.string());

      fout_nis << header << "\n\n";

      fout_nis.precision(12);
      fout_nis.setf(ios::fixed, ios::floatfield);

      fout_nis << "E_lab,E_cm" << endl;

      list<double> ni_energies =
          bqptr->getFreeTwoParticleEnergiesInElab(emin, emax);
      for (double& ni_energy : ni_energies) {
        double ecm_energy = bqptr->getEcmOverMrefFromElab(ni_energy);
        fout_nis << ni_energy << "," << ecm_energy << endl;
      }
      fout_nis.close();

      double einc_percent = elabs_inc[blocknum];
      double elab = emin;

      ni_energies.push_front(emin);
      ni_energies.push_back(emax);

      vector<double> elabvals;
      for (int i = 0; i < ni_energies.size() - 1; ++i) {
        double e1 = *std::next(ni_energies.begin(), i);
        double e2 = *std::next(ni_energies.begin(), i + 1);
        double einc = einc_percent * (e2 - e1); // relative to NI interval
        while (elab <= e2 + 1e-9) {
          elabvals.push_back(elab);
          elab += einc;
        }
      }

      uint nvals = elabvals.size();
      if (outmode == "full")
        nsamp = 0;
      vector<CVector> omegavals(nvals, CVector(nsamp + 1));
      vector<vector<CVector>> eigenvals(
          bqptr->getBasisSize(), vector<CVector>(nvals, CVector(nsamp + 1)));

      for (uint b = 0; b <= nsamp; ++b) {
        CMatrix last_iter_eigenvectors;
        for (uint k = 0; k < nvals; ++k) {
          omegavals[k][b] =
              bqptr->getOmegaFromElab(omega_mu, elabvals[k], qctype_enum);
          if (do_print_eigenvals) {
            CVector ev_res;
            ev_res = bqptr->getQCEigenvaluesFromElab(elabvals[k], qctype_enum);
            for (int dim = 0; dim < bqptr->getBasisSize(); ++dim) {
              eigenvals[dim][k][b] = ev_res[dim];
            }
          }
        }
      }

      vector<double> ecmvals;
      for (uint k = 0; k < nvals; ++k) {
        double ecm_energy = bqptr->getEcmOverMrefFromElab(elabvals[k]);
        ecmvals.push_back(ecm_energy);
      }

      ofstream fout_omega(omega_filepath.string());

      fout_omega << header << "\n\n";

      fout_omega.precision(12);
      fout_omega.setf(ios::fixed, ios::floatfield);
      fout_omega << "E_lab,E_cm,Omega_re,Omega_im" << endl;
      if (nsamp == 0) {
        for (uint k = 0; k < nvals; ++k)
          fout_omega << elabvals[k] << "," << ecmvals[k] << ","
                     << omegavals[k][0].real() << "," << omegavals[k][0].imag()
                     << endl;
      }

      if (do_root_find) {
        std::vector<double> roots;
        std::vector<uint> fn_calls;
        bqptr->getEcmRootsInElabInterval(omega_mu, emin, emax, qctype_enum,
                                         root_config, roots, fn_calls);
        ofstream fout_roots(roots_filepath.string());

        fout_roots << header << "\n\n";

        fout_roots.precision(12);
        fout_roots.setf(ios::fixed, ios::floatfield);

        fout_roots << "E_lab,E_cm" << endl;

        if (nsamp == 0) {
          for (double root : roots) {
            double root_ecm = bqptr->getEcmOverMrefFromElab(root);
            fout_roots << root << "," << root_ecm << endl;
          }
        }

        // put fn calls in the log
        logger << "Number of function calls for each root finding interval:"
               << endl;
        for (int i = 0; i < fn_calls.size(); ++i) {
          logger << "Interval " << i << ": " << fn_calls[i] << endl;
        }
      }
      // } else if (m_obs->isJackknifeMode()) {
      //   MCEstimate mcest;
      //   for (uint k = 0; k < nvals; ++k) {
      //     m_obs->jack_analyze(omegavals[k], mcest);
      //     fout_omega << "  " << setw(12) << elabvals[k] << " " << setw(20)
      //                << mcest.getAverageEstimate() << " " << setw(20)
      //                << mcest.getSymmetricError() << endl;
      //   }
      // } else {
      //   MCEstimate mcest;
      //   for (uint k = 0; k < nvals; ++k) {
      //     m_obs->boot_analyze(omegavals[k], mcest);
      //     double avg = mcest.getAverageEstimate();
      //     double upperr = mcest.getUpperConfLimit() - avg;
      //     double dwnerr = mcest.getLowerConfLimit() - avg;
      //     fout_omega << "  " << setw(12) << elabvals[k] << " " << setw(20)
      //                << avg << " " << setw(20) << upperr << " " << setw(20)
      //                << dwnerr << endl;
      //   }
      // } // need to add CVector support for jackknife and bootstrap.
      // Probably just look at SigMond
      fout_omega.close();

      if (do_print_eigenvals) {
        ofstream fout_eigenvals(eigenvals_filepath.string());
        fout_eigenvals << header << "\n\n";

        fout_eigenvals.precision(12);
        fout_eigenvals.setf(ios::fixed, ios::floatfield);

        fout_eigenvals << "E_lab,E_cm";
        for (int dim = 0; dim < bqptr->getBasisSize(); ++dim) {

          fout_eigenvals << ",ev_re,ev_im" << dim;
          // if (m_obs->isJackknifeMode() && nsamp != 0) {
          //   fout_eigenvals << "_AverageEstimate,ev" << dim <<
          //   "_SymmetricError";
          // } else if (m_obs->isBootstrapMode() && nsamp != 0) {
          //   fout_eigenvals << "_AverageEstimate,ev" << dim <<
          //   "_UpperError,ev"
          //                  << dim << "_LowerError";
          // }
        }
        fout_eigenvals << endl;

        if (nsamp == 0) {
          for (uint k = 0; k < nvals; ++k) {
            fout_eigenvals << elabvals[k] << "," << ecmvals[k];
            for (int dim = 0; dim < bqptr->getBasisSize(); ++dim) {
              fout_eigenvals << "," << eigenvals[dim][k][0].real() << ","
                             << eigenvals[dim][k][0].imag();
            }
            fout_eigenvals << endl;
          }
        }
        // } else if (m_obs->isJackknifeMode()) {
        //   MCEstimate mcest;
        //   for (uint k = 0; k < nvals; ++k) {
        //     for (int dim = 0; dim < bqptr->getBasisSize(); ++dim) {
        //       m_obs->jack_analyze(eigenvals[dim][k], mcest);
        //       fout_omega << "," << mcest.getAverageEstimate() << ","
        //                  << mcest.getSymmetricError();
        //     }
        //     fout_omega << endl;
        //   }
        // } else {
        //   MCEstimate mcest;
        //   for (uint k = 0; k < nvals; ++k) {
        //     for (int dim = 0; dim < bqptr->getBasisSize(); ++dim) {
        //       m_obs->boot_analyze(eigenvals[dim][k], mcest);
        //       double avg = mcest.getAverageEstimate();
        //       double upperr = mcest.getUpperConfLimit() - avg;
        //       double dwnerr = mcest.getLowerConfLimit() - avg;
        //       fout_omega << avg << "," << upperr << "," << dwnerr;
        //     }
        //     fout_omega << endl;
        //   }
        // }
        fout_eigenvals.close();
      }
    }

    m_obs->clearSamplings();
    XMLHandler xmlK;
    if (Kmat != 0)
      Kmat->output(xmlK);
    else
      Kinv->output(xmlK);
    xmlout.put_child(xmlK);
    XMLHandler xmllog;
    xmlformat("Logger", logger.str(), xmllog);
    xmlout.put_child(xmllog);
  } catch (const std::exception& xp) {
    string msg("doPrint failed: ");
    msg += xp.what();
    cout << msg << endl;
    throw(std::runtime_error(msg));
  }
}

// ***************************************************************************************
