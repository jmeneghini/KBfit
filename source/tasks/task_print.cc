#include "chisq_detres.h"
#include "chisq_fit.h"
#include "task_handler.h"
#include "root_finder.h"
#include <iomanip>

using namespace std;

// **************************************************************************
// *                                                                        *
// *    Prints out the Omega function for Elab between a minimum and        *
// *    maximum specified by user at fixed increment, specified by user.    *
// *    Output for each KBBlock is sent to a separate file with a common    *
// *    stub and different integer suffices.                                *
// *    Format of output files:                                             *
// *     -- if OutputMode = "full", format is                               *
// *                     Elab/mref  value_from_full_sample                  *
// *     -- if OutputMode = "resampled", format is                          *
// *              Elab/mref  value_from_full_sample upper_error down_error  *
// *                                                                        *
// *    Notes:                                                              *
// *                                                                        *
// *    -- An L^3 spatial lattice is required.  Eventually, the length      *
// *       times the reference scale is needed.  To determine this,         *
// *       the reference scale times the temporal lattice spacing must      *
// *       be specified in <ReferenceMassTimeSpacingProduct>.               *
// *                                                                        *
// *    -- If <OmegaMu> is specified, then the Omega function is used       *
// *       in the residual.  If absent, the determinant itself is used.     *
// *                                                                        *
// *    -- Either <KtildeMatrixInverse> or <KtildeMatrix> must be given.    *
// *       Depending on which is input, det(1-K*B) or det(K^(-1)-B)         *
// *       is used.                                                         *
// *                                                                        *
// *    -- If using an anisotropic lattice, a tag <LatticeAnisotropy>,      *
// *       which is the spatial over the temporal spacing, must be given.   *
// *       If this tag is absent, an isotropic lattice is assumed.          *
// *                                                                        *
// *    -- All lab-frame energies and particle masses can be input either   *
// *       as ratios of a reference mass or as a product with the time      *
// *       spacing of the lattice.  The default format should be specified  *
// *       in the tag <DefaultEnergyFormat> whose value can be either       *
// *       "reference_ratio" or "time_spacing_product".                     *
// *                                                                        *
// *    Format of input XML:                                                *
// *                                                                        *
// *                                                                        *
// *    <Task>                                                              *
// *                                                                        *
// *     <Action>DoPrint</Action>                                           *
// *                                                                        *
// *      <OutputStub>...</OutputStub>                                      *
// *                                                                        *
// *      <OutputMode>full</OutputMode> or "resampled"                      *
// *                                                                        *
// *      <OmegaMu>8.0</OmegaMu>  (optional)                                *
// *                                                                        *
// *      <QuantizationCondition>...</QuantizationCondition>                *
// *       (StildeCB, StildeinvCB, KtildeB, KtildeinvB)                     *
// *                                                                        *
// *  <RootFinder>                                                          *
// *    <LabFrameEnergyMin>1.10</LabFrameEnergyMin>                         *
// *    <LabFrameEnergyMax>2.30</LabFrameEnergyMax>                         *
// *                                                                        *
// *    <AdaptiveBracket>                                                   *
// *                                                                        *
// *      <!--‑‑‑‑‑ All tags below are OPTIONAL ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑>*
// *                                                                        *
// *      <InitialStepSize>         1e-2  </InitialStepSize>                *
// *      <AbsXTolerance>           1e-6  </AbsXTolerance>                  *
// *      <AbsResidualTolerance>    1e-10 </AbsResidualTolerance>           *
// *      <MinStepSize>             1e-4  </MinStepSize>                    *
// *      <MaxStepSize>             5e-2  </MaxStepSize>                    *
// *      <StepScaleLimit>          3.0   </StepScaleLimit>                 *
// *      <PlateauMod2Threshold>    1e-4  </PlateauMod2Threshold>           *
// *      <PlateauCountBeforeJump>  4     </PlateauCountBeforeJump>         *
// *                                                                        *
// *    </AdaptiveBracket>                                                  *
// *  </RootFinder>                                                         *
// *                                                                        *
// *----------------------------------------------------------------------  *
// *            TAG  →  AdaptiveBracketConfig FIELD MAPPING                 *
// *----------------------------------------------------------------------  *
// *      <InitialStepSize>          →  initial_step_size                   *
// *      <AbsXTolerance>            →  x_tol                               *
// *      <AbsResidualTolerance>     →  zero_tol                            *
// *      <MinStepSize>              →  min_step_size                       *
// *      <MaxStepSize>              →  max_step_size                       *
// *      <StepScaleLimit>           →  step_scale_limit                    *
// *      <PlateauMod2Threshold>     →  plateau_mod2_threshold              *
// *      <PlateauCountBeforeJump>   →  plateau_count_before_jump           *
// *------------------------------------------------------------------------*
// *                                                                        *
// *      <KtildeMatrixInverse>  (or <KtildeMatrix> - match quant cond)     *
// *          .....                                                         *
// *      </KtildeMatrixInverse>                                            *
// *                                                                        *
// *      <DefaultEnergyFormat>reference_ratio</DefaultEnergyFormat>        *
// *                   (or time_spacing_product)                            *
// *                                                                        *
// *      <MCEnsembleParameters>...</MCEnsembleParameters>                  *
// *        ... one for each Monte Carlo ensemble                           *
// *                                                                        *
// *      <KBBlock>...</KBBlock>                                            *
// *        ... one for each KB quantization block                          *
// *                                                                        *
// *      <KBObservables>                                                   *
// *       <MCSamplingInfo>...</MCSamplingInfo>                             *
// *       <SamplingData>                                                   *
// *          <FileName>...</FileName>  (all sampling files needed to       *
// *             ....                    obtain all above MCObsInfo's)      *
// *       </SamplingData>                                                  *
// *      </KBObservables>                                                  *
// *                                                                        *
// *    </Task>                                                             *
// *                                                                        *
// *                                                                        *
// *                                                                        *
// *    For each Monte Carlo ensemble involved, there should be a tag       *
// *    with the format:                                                    *
// *                                                                        *
// *      <MCEnsembleParameters>                                            *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <ReferenceMassTimeSpacingProduct>                               *
// *            <MCObs>...</MCObs>                                          *
// *        </ReferenceMassTimeSpacingProduct>                              *
// *        <LatticeAnisotropy>        (optional a_s/a_t) (unity if absent) *
// *            <MCObs>...</MCObs>                                          *
// *        </LatticeAnisotropy>                                            *
// *        <ParticleMass>                                                  *
// *           <Name>pion</Name> (should match names used in K-matrix)      *
// *           <MCObs>...</MCObs>                                           *
// *        </ParticleMass>                                                 *
// *           ... other particle masses                                    *
// *      </MCEnsembleParameters>                                           *
// *                                                                        *
// *    When specifying an MCObsInfo, either a short form or a long form    *
// *    can be used (must be nonsimple and real):                           *
// *                                                                        *
// *     <MCObservable>                                                     *
// *       <ObsName>T1up_Energy</ObsName> (32 char or less, no blanks)      *
// *       <Index>3</Index>        (opt nonneg integer: default 0)          *
// *     </MCObservable>                                                    *
// *                                                                        *
// *     <MCObs>T1up_Energy 3</MCObs>                                       *
// *                                                                        *
// *                                                                        *
// *    For each KB quantization block, a tag of the form is needed:        *
// *                                                                        *
// *      <KBBlock>                                                         *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <BoxQuantization>                                               *
// *          <TotalMomentumRay>ar</TotalMomentumRay>                       *
// *          <TotalMomentumIntSquared>0</TotalMomentumIntSquared>          *
// *          <LGIrrep>T1u</LGIrrep>                                        *
// *          <LmaxValues>5 3</LmaxValues>  (one for each decay channel)    *
// *        </BoxQuantization>                                              *
// *        <LabFrameEnergyMin>...</LabFrameEnergyMin>                      *
// *        <LabFrameEnergyMax>...</LabFrameEnergyMax>                      *
// *        <LabFrameEnergyInc>...</LabFrameEnergyInc>                      *
// *      </KBBlock>                                                        *
// *                                                                        *
// *                                                                        *
// *                                                                        *
// *    Specification of the K-matrix or the inverse of the K-matrix is     *
// *    done using XML of the format:                                       *
// *                                                                        *
// *      <KtildeMatrixInverse>  (or <KtildeMatrix>)                        *
// *                                                                        *
// *        <DecayChannels>                                                 *
// *           <DecayChannelInfo>                                           *
// *              <Particle1Name>pion</Particle1Name>                       *
// *              <Spin1TimesTwo>0</Spin1TimesTwo>                          *
// *              <Identical/> (if identical, do not include tags below)    *
// *              <Particle2Name>eta</Particle2Name>                        *
// *              <Spin2TimesTwo>2</Spin2TimesTwo>                          *
// *             <IntrinsicParities>same</IntrinsicParities> (or "opposite")*
// *           </DecayChannelInfo>                                          *
// *                                                                        *
// *            ... other channels infos ...                                *
// *                                                                        *
// *          (Order matters: first <DecayChannelInfo> tag is channel 0,    *
// *           second <DecayChannelInfo> is channel 1, and so on.  In the   *
// *           K-matrix, channels are referred to using the index 0, 1,...) *
// *                                                                        *
// *        </DecayChannels>                                                *
// *                                                                        *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>2</JTimesTwo>                                    *
// *            <KIndex>L(1) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(1) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *             <Polynomial><Powers>1 3</Powers></Polynomial>              *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>6</JTimesTwo>                                    *
// *            <KIndex>L(3) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(3) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *             <Polynomial><Degree>0 </Degree></Polynomial>               *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *        <Element>                                                       *
// *          <KElementInfo>                                                *
// *            <JTimesTwo>10</JTimesTwo>                                   *
// *            <KIndex>L(5) 2S(0) chan(0)</KIndex>                         *
// *            <KIndex>L(5) 2S(0) chan(0)</KIndex>                         *
// *          </KElementInfo>                                               *
// *          <FitForm>                                                     *
// *            <Polynomial><Degree>0</Degree></Polynomial>                 *
// *          </FitForm>                                                    *
// *        </Element>                                                      *
// *                                                                        *
// *        <StartingValues>                                                *
// *                                                                        *
// *          <KFitParamInfo>                                               *
// *            <PolynomialTerm>                                            *
// *               <Power>3</Power>                                         *
// *               <KElementInfo>...</KElementInfo>                         *
// *            </PolynomialTerm>                                           *
// *            <StartingValue>0.4534</StartingValue>                       *
// *          </KFitParamInfo>                                              *
// *          <KFitParamInfo>                                               *
// *            <PoleEnergy>                                                *
// *               <Index>3</Index>                                         *
// *               <JTimesTwo>2</JTimesTwo>                                 *
// *            </PoleEnergy>                                               *
// *            <StartingValue>2.2</StartingValue>                          *
// *          </KFitParamInfo>                                              *
// *          <KFitParamInfo>                                               *
// *            <PoleCoupling>                                              *
// *               <Index>3</Index>                                         *
// *               <JTimesTwo>2</JTimesTwo>                                 *
// *               <KIndex>...</KIndex>                                     *
// *            </PoleCoupling>                                             *
// *            <StartingValue>3.3</StartingValue>                          *
// *          </KFitParamInfo>                                              *
// *                                                                        *
// *               ....                                                     *
// *        </StartingValues>                                               *
// *                                                                        *
// *      </KtildeMatrixInverse>                                            *
// *                                                                        *
// *                                                                        *
// *    Final note: The reference mass, the particle masses, and the        *
// *    anisotropy for each ensemble can be set to a fixed value by         *
// *    replacing the <MCObs>/<MCObservable> tag by a                       *
// *      <FixedValue>1.1123</FixedValue> tag.                              *
// *    Some groups, such as JLab, erroneously use such fixed values.       *
// *    The above feature is useful for determining how such an erroneous   *
// *    procedure effects the final error estimates on the K-matrix fix     *
// *    parameters. (NOTE: the fixed values for particle masses always      *
// *    refer to energy ratios over the reference mass.)                    *
// *                                                                        *
// *                                                                        *
// **************************************************************************

void TaskHandler::doPrint(XMLHandler& xmltask, XMLHandler& xmlout,
                          int taskcount) {
  try {
    xmlout.set_root("PrintOmega");
    stringstream logger;
    logger.precision(12);
    KtildeMatrixCalculator* Kmat = 0;
    KtildeInverseCalculator* Kinv = 0;

    string outstub;
    xmlread(xmltask, "OutputStub", outstub, "doPrint");
    outstub = tidyString(outstub);

    // check if we're printing eigenvalues
    int count_print_eigenvals =
        xmltask.count_among_children("PrintEigenvalues");
    if (count_print_eigenvals > 1) {
      throw(std::invalid_argument(
          "Multiple PrintEigenvalues tags cannot be present"));
    }
    bool do_print_eigenvals = count_print_eigenvals;

    if (outstub.empty())
      throw(std::runtime_error("No output stub specified"));
    string outmode = "full";
    xmlreadif(xmltask, "OutputMode", outmode, "doPrint");
    if ((outmode != "full") && (outmode != "resampled"))
      throw(std::runtime_error("Invalid output mode specified"));

    //  get sampling mode and names of input files
    XMLHandler xmlr(xmltask, "KBObservables");
    MCSamplingInfo sampinfo(xmlr);
    if (sampinfo != m_obs->getSamplingInfo()) {
      // m_obs->setSamplingInfo(sampinfo);
      // logger << "MCSamplingInfo reset in KBObsHandler in
      // DeterminantResidualFit"<<endl;}
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
    xmlreadif(xmltask, "OmegaMu", omega_mu, "DeterminantResidualFit");
    logger << "Omega mu = " << omega_mu << endl;

    //  get format for energies/mass (in terms of product with
    //  lattice time spacing or a ratio with reference mass)
    bool energy_ratios = true;
    {
      string reply;
      xmlreadif(xmltask, "DefaultEnergyFormat", reply,
                "DeterminantResidualFit");
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
      throw(std::invalid_argument(
              "Multiple RootFinder tags cannot be present"));
    bool do_root_find = root_counts;
    AdaptiveBracketConfig root_config;
    if (do_root_find) {
      XMLHandler xmlroot(xmltask, "RootFinder");
      root_config = AdaptiveBracketRootFinder::makeConfigFromXML(xmlroot);
    }


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
      DeterminantResidualFit::read_obs(*it, "ReferenceMassTimeSpacingProduct",
                                       false, rkey, kset, pname, mcens,
                                       fixed_values);
      ref_at_mass.insert(make_pair(mcens, rkey));
      if (!energy_ratios) {
        kset.erase(rkey);
      } // remove ref_at_mass from kset
      // get anisotropy info
      if (xml_tag_count(*it, "LatticeAnisotropy") == 1) {
        DeterminantResidualFit::read_obs(*it, "LatticeAnisotropy", false, rkey,
                                         kset, pname, mcens, fixed_values);
        anisotropy.insert(make_pair(mcens, rkey));
      }
      // get particle mass infos
      map<string, MCObsInfo> pmap;
      list<XMLHandler> xmlp = it->find("ParticleMass");
      for (list<XMLHandler>::iterator pt = xmlp.begin(); pt != xmlp.end();
           ++pt) {
        DeterminantResidualFit::read_obs(*pt, "ParticleMass", true, rkey, kset,
                                         pname, mcens, fixed_values);
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
      if (qctype_enum == BoxQuantization::StildeCB ||
          qctype_enum == BoxQuantization::KtildeB) {
        if (k2 == 1) {
          throw(std::invalid_argument(
              "KtildeMatrixInverse cannot be used with StildeCB or KtildeB"));
        }
      }
      if (qctype_enum == BoxQuantization::StildeinvCB ||
          qctype_enum == BoxQuantization::KtildeinvB) {
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

    //  Create a folder with project name if it does not already exist,
    //  and inside that folder creat a folder with the QuantCond,
    //  where the output files will be stored.
    string project_name = outstub;
    string quant_cond_name = qctype;

    filesystem::path project_dir = filesystem::path(project_name);
    filesystem::path quant_cond_dir = project_dir / quant_cond_name;

    std::error_code ec;
    if (!filesystem::exists(quant_cond_dir)) {
      if (!filesystem::create_directories(quant_cond_dir, ec)) {
        if (ec) {
          throw(
              std::runtime_error("Error creating directory: " + ec.message()));
        }
      }
    }
    // Now move into folder
    filesystem::current_path(quant_cond_dir);

    //  Now start the output, block by block

    for (uint blocknum = 0; blocknum < BQ.size(); ++blocknum) {
      const MCEnsembleInfo& mcens = blockens[blocknum];
      string omega_filename = "Omega_Block" + make_string(blocknum) + ".csv";
      string eigenvals_filename =
          "Eigenvals_Block" + make_string(blocknum) + ".csv";
      string nis_filename =
          "NonInteractingEnergies_Block" + make_string(blocknum) + ".csv";
      string roots_filename = "Roots_Block" + make_string(blocknum) + ".csv";

      BoxQuantization* bqptr = BQ[blocknum];
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

      double emin = elabs_min[blocknum];
      double emax = elabs_max[blocknum];
      double einc = elabs_inc[blocknum];
      double elab = emin;
      vector<double> elabvals;
      while (elab <= emax) {
        elabvals.push_back(elab);
        elab += einc;
      }
      uint nvals = elabvals.size();
      if (outmode == "full")
        nsamp = 0;
      vector<CVector> omegavals(nvals, CVector(nsamp + 1));
      vector<vector<CVector>> eigenvals(
          bqptr->getBasisSize(), vector<CVector>(nvals, CVector(nsamp + 1)));

      for (uint b = 0; b <= nsamp; ++b) {
        bqptr->setRefMassL(mrefL[b]);
        for (uint ci = 0; ci < bqptr->getNumberOfDecayChannels(); ++ci) {
          bqptr->setMassesOverRef(ci, (*(particlemass1[ci]))[b],
                                  (*(particlemass2[ci]))[b]);
        }
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

      ofstream fout_omega(omega_filename);

      string header = "#" + mcens.str() + " # MomRay " + bqptr->getMomRay() +
                      " # P^2 = " +
                      std::to_string(bqptr->getTotalMomentumIntegerSquared()) +
                      " # Box Irrep " + bqptr->getLittleGroupBoxIrrep();

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
        bqptr->getRootsInElabInterval(omega_mu, emin, emax, qctype_enum,
                                      root_config, roots);
        ofstream fout_roots(roots_filename);
        fout_roots << header << "\n\n";

        fout_roots.precision(12);
        fout_roots.setf(ios::fixed, ios::floatfield);

        fout_roots << "E_lab,E_cm";

        if (nsamp == 0) {
          for (uint k = 0; k < nvals; ++k) {
            double root_ecm = bqptr->getEcmOverMrefFromElab(roots[k]);
            fout_roots << roots[k] << "," << root_ecm << endl;
          }
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
        ofstream fout_eigenvals(eigenvals_filename);
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
        fout_omega.close();
      }
      ofstream fout_nis(nis_filename);

      fout_nis << header << "\n\n";

      fout_nis.precision(12);
      fout_nis.setf(ios::fixed, ios::floatfield);

      fout_nis << "E_lab,E_cm" << endl;

      list<double> ni_energies = bqptr->getFreeTwoParticleEnergies(emin, emax);
      for (double& ni_energy : ni_energies) {
        double ecm_energy = bqptr->getEcmOverMrefFromElab(ni_energy);
        fout_nis << ni_energy << "," << ecm_energy << endl;
      }
      fout_nis.close();
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
