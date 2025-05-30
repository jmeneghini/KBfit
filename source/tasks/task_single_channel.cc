#include "K_matrix_info.h"
#include "box_matrix.h"
#include "chisq_detres.h"
#include "chisq_fit.h"
#include "sampling_info.h"
#include "task_handler.h"
#include "xml_handler.h"
#include <filesystem>
#include <iomanip>
#include <string>

using namespace std;

// **************************************************************************
// *     Prints out q_a/mref cot(delta_0), where a is the specified decay   *
// *     channel and delta_0 is the phase shift of the s-wave.              *
// *     The calculation holds for elastic scattering between two spinless  *
// *     particles. Additionally, the s-wave must subduce onto the provided *
// *     box irreps.                                                        *
// *                                                                        *
// *     Format of output files:                                            *
// *     -- if OutputMode = "full", format is                               *
// *              Elab/mref  value_from_full_sample                         *
// *     -- if OutputMode = "resampled", format is                          *
// *              Elab/mref  value_from_full_sample upper_error down_error  *
// *                                                                        *
// *   Notes:                                                               *
// *    -- An L^3 spatial lattice is required.  Eventually, the length      *
// *       times the reference scale is needed.  To determine this,         *
// *       the reference scale times the temporal lattice spacing must      *
// *       be specified in <ReferenceMassTimeSpacingProduct>.               *
// *                                                                        *
// *    -- All lab-frame energies and particle masses can be input either   *
// *       as ratios of a reference mass or as a product with the time      *
// *       spacing of the lattice.  The default format should be specified  *
// *       in the tag <DefaultEnergyFormat> whose value can be either       *
// *       "reference_ratio" or "time_spacing_product".                     *
// *                                                                        *
// *    Format of input XML:                                                *
// *                                                                        *
// *    <Task>                                                              *
// *     <Action>SingleChannel</Action>                                     *
// *                                                                        *
// *      <OutputStub>...</OutputStub>                                      *
// *                                                                        *
// *      <OutputMode>full</OutputMode> or "resampled"                      *
// *                                                                        *
// *      <DefaultEnergyFormat>reference_ratio</DefaultEnergyFormat>        *
// *                   (or time_spacing_product)                            *
// *                                                                        *
// *      <MCEnsembleParameters>...</MCEnsembleParameters>                  *
// *        ... one for each Monte Carlo ensemble                           *
// *                                                                        *
// *      <KBBlock>...</KBBlock>                                            *
// *        ... one for each box irrep                                      *
// *                                                                        *
// *      <DecayChannels>                                                   *
// *        <DecayChannelInfo>                                              *
// *            <Particle1Name>pion</Particle1Name>                         *
// *            <Spin1TimesTwo>0</Spin1TimesTwo>                            *
// *            <Identical/> (if identical, do not include tags below)      *
// *            <Particle2Name>kaon</Particle2Name>                         *
// *            <Spin2TimesTwo>0</Spin2TimesTwo>                            *
// *            <IntrinsicParities>same</IntrinsicParities> (or "opposite") *
// *        </DecayChannelInfo>                                             *
// *      </DecayChannels>                                                  *
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
// *    For each Monte Carlo ensemble involved, there should be a tag       *
// *    with the format:                                                    *
// *                                                                        *
// *      <MCEnsembleParameters>                                            *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <ReferenceMassTimeSpacingProduct>                               *
// *            <MCObs>...</MCObs>                                          *
// *        </ReferenceMassTimeSpacingProduct>                              *
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
// *
// *                                                                        *
// *      <KBBlock>                                                         *
// *        <MCEnsembleInfo>...</MCEnsembleInfo>                            *
// *        <BoxQuantization>                                               *
// *          <TotalMomentumRay>ar</TotalMomentumRay>                       *
// *          <TotalMomentumIntSquared>0</TotalMomentumIntSquared>          *
// *          <LGIrrep>A1g</LGIrrep>                                        *
// *          <LmaxValues>0</LmaxValues>      <-- exact requirement         *
// *        </BoxQuantization>                                              *
// *        <LabFrameEnergy>                                                *
// *            <MCObs>...</MCObs>                                          *
// *        </LabFrameEnergy>                                               *
// *      </KBBlock>                                                        *
// *                                                                        *

void TaskHandler::doSingleChannel(XMLHandler& xmltask, XMLHandler& xmlout,
                                  int taskcount) {
  try {
    xmlout.set_root("SingleChannel");
    stringstream logger;
    logger.precision(12);

    string outstub;
    xmlread(xmltask, "OutputStub", outstub, "doSingleChannel");
    outstub = tidyString(outstub);

    if (outstub.empty())
      throw(std::runtime_error("No output stub specified"));
    string outmode = "full";
    xmlreadif(xmltask, "OutputMode", outmode, "doSingleChannel");
    if ((outmode != "full") && (outmode != "resampled"))
      throw(std::runtime_error("Invalid output mode specified"));

    //  get sampling mode and names of input files
    XMLHandler xmlr(xmltask, "KBObservables");
    MCSamplingInfo sampinfo(xmlr);
    if (sampinfo != m_obs->getSamplingInfo()) {
      throw(std::invalid_argument("KBObservables MCSamplingInfo does not match "
                                  "that of the KBObsHandler"));
    }

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

    // get the decay channel info
    XMLHandler xmld(xmltask, "DecayChannels");
    list<XMLHandler> dcxml = xmld.find_among_children("DecayChannelInfo");
    if (dcxml.size() != 1) {
      throw(std::invalid_argument(
          "SingleChannel task requires exactly one DecayChannelInfo tag"));
    }

    // ensure that the particles are spinless and store the channel info
    DecayChannelInfo single_channel(dcxml.front());
    if (single_channel.getSpin1timestwo() != 0) {
      throw(std::invalid_argument(
          "Spin1TimesTwo must be 0 for spinless particles"));
    }
    if (single_channel.getSpin2timestwo() != 0) {
      throw(std::invalid_argument(
          "Spin2TimesTwo must be 0 for spinless particles"));
    }

    // get set of particle names from decay channel
    set<string> pnames;
    pnames.insert(single_channel.getParticle1Name());
    pnames.insert(single_channel.getParticle2Name());

    // get the energy format (the independent variable associated with each
    // cot(delta_0))
    bool energy_ratios = true;
    {
      string reply;
      xmlreadif(xmltask, "DefaultEnergyFormat", reply, "doSingleChannel");
      if ((reply == "reference_ratio") || (reply.empty()))
        energy_ratios = true;
      else if (reply == "time_spacing_product")
        energy_ratios = false;
      else
        throw(std::invalid_argument("Invalid <DefaultEnergyFormat> tag"));
    }

    // Loop over the different MC ensembles to get information
    // about all needed parameters.  Particle masses, anisotropy,
    // parameters should be associated with the ensemble.
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

    //  Loop over the KB quantization blocks to get the lab-frame energies
    //  infos.
    vector<KBObsInfo> labenergies;
    vector<MCEnsembleInfo> blockens;
    map<MCEnsembleInfo, uint> ensemble_idmap;
    vector<BoxQuantization*> BQ;
    vector<vector<uint>>
        block_energy_indices; // Track which energies belong to which block

    list<XMLHandler> xmlkb = xmltask.find("KBBlock");
    uint blockcount = 0;
    uint ensemblecount = 0;

    // dummy KtildeMatrixCalculator pointer for BoxQuantization
    KElementInfo single_channel_info(0, 0, 0, 0, 0, 0, 0);
    std::list<std::pair<KElementInfo, Polynomial>> dummy_pelems = {
        pair<KElementInfo, Polynomial>(single_channel_info, Polynomial())};
    std::list<std::pair<KElementInfo, SumOfPoles>> dummy_selems = {};
    std::list<std::pair<KElementInfo, SumOfPolesPlusPolynomial>> dummy_spelems =
        {};
    KtildeMatrixCalculator* dummyKmat = new KtildeMatrixCalculator(
        dummy_pelems, dummy_selems, dummy_spelems, {single_channel});

    for (list<XMLHandler>::iterator it = xmlkb.begin(); it != xmlkb.end();
         ++it) {
      // construct the BoxQuantization object for this block
      BoxQuantization* bqptr = new BoxQuantization(*it, dummyKmat, nullptr);
      BQ.push_back(bqptr);

      MCEnsembleInfo mcens(*it);
      if (ensembles.find(mcens) == ensembles.end())
        throw(
            std::invalid_argument("KBBlock associated with unknown ensemble"));
      blockens.push_back(mcens);
      if (ensemble_idmap.find(mcens) == ensemble_idmap.end()) {
        ensemble_idmap.insert(make_pair(mcens, ensemblecount));
        ensemblecount++;
      }
      set<MCObsInfo>& kset = needed_keys[mcens];
      uint nres = 0;
      MCObsInfo rkey;
      vector<uint>
          this_block_energies; // Track energies for this specific block
      list<XMLHandler> xmlee = it->find("LabFrameEnergy");
      for (list<XMLHandler>::iterator eet = xmlee.begin(); eet != xmlee.end();
           ++eet) {
        DeterminantResidualFit::read_obs(*eet, "LabFrameEnergy", rkey, kset);
        uint energy_index = labenergies.size();
        labenergies.push_back(KBObsInfo(mcens, rkey));
        this_block_energies.push_back(energy_index);
        nres++;
      }
      if (nres == 0)
        throw(std::invalid_argument(
            "No energies available in at least one block"));
      block_energy_indices.push_back(this_block_energies);
      blockcount++;
    }
    if (blockcount == 0) {
      throw(std::runtime_error("No data to analyze"));
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
          const RVector& anisotropy_values =
              m_obs->getFullAndSamplingValues(anisotropykey);
          if (anisotropy_values.size() != (nsamp + 1))
            throw(std::runtime_error("Resampling size mismatch in KBfit"));
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat) * anisotropy_values[k];
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

    logger << "Found " << labenergies.size() << " energy observables across "
           << blockcount << " blocks" << endl;

    //  Create a folder with project name if it does not already exist.
    //  The output files will be stored here.
    string project_name = outstub;

    filesystem::path project_dir = filesystem::path(project_name);

    std::error_code ec;
    if (!filesystem::exists(project_dir)) {
      if (!filesystem::create_directories(project_dir, ec)) {
        if (ec) {
          throw(
              std::runtime_error("Error creating directory: " + ec.message()));
        }
      }
    }
    // Now move into folder
    filesystem::current_path(project_dir);

    // Now start the single-channel computations, block by block (following
    // task_print.cc pattern)
    for (uint blocknum = 0; blocknum < BQ.size(); ++blocknum) {
      const MCEnsembleInfo& mcens = blockens[blocknum];

      string filename = "Block" + make_string(blocknum) + ".csv";

      BoxQuantization* bqptr = BQ[blocknum];

      logger << "Filename = " << filename << endl;
      logger << "Block " << blocknum << ": " << mcens.str() << endl;
      logger << "MomRay " << bqptr->getMomRay()
             << "   P^2 = " << bqptr->getTotalMomentumIntegerSquared()
             << " Box Irrep " << bqptr->getLittleGroupBoxIrrep() << endl;
      uint nsamp = m_obs->getNumberOfResamplings();
      map<string, MCObsInfo>& pmap = particle_masses[mcens];

      KBObsInfo lengthkey(mcens, MCObsInfo("LengthReference"));
      const RVector& mrefL = m_obs->getFullAndSamplingValues(lengthkey);

      // Single channel logic - use the channel info from the beginning
      const string& pname1 = single_channel.getParticle1Name();
      const string& pname2 = single_channel.getParticle2Name();
      KBObsInfo mass1key(mcens, pmap[pname1]);
      KBObsInfo mass2key(mcens, pmap[pname2]);
      const RVector* particlemass1 =
          &(m_obs->getFullAndSamplingValues(mass1key));
      const RVector* particlemass2 =
          &(m_obs->getFullAndSamplingValues(mass2key));

      // Set masses for all resamplings - single channel only
      for (int b = 0; b <= nsamp; ++b) {
        bqptr->setRefMassL(mrefL[b]);
        bqptr->setMassesOverRef(0, (*particlemass1)[b], (*particlemass2)[b]);
      }

      string header = "#" + mcens.str() + " # MomRay " + bqptr->getMomRay() +
                      " # P^2 = " +
                      std::to_string(bqptr->getTotalMomentumIntegerSquared()) +
                      " # Box Irrep " + bqptr->getLittleGroupBoxIrrep();

      ofstream fout(filename);

      fout << header << "\n\n";

      fout.precision(12);
      fout.setf(ios::fixed, ios::floatfield);

      // Give results for E_lab/mref, E_cm/mref, and (q/mref)^2

      if (outmode == "full") {
        fout << "E_lab/mref,E_cm/mref,(q/mref)^2,q/mref_cot_delta" << endl;
      } else {
        fout << "E_lab/mref,E_lab/mref_upper_err,E_lab/mref_lower_err,E_lab/"
                "mref_sym_err,"
             << "E_cm/mref,E_cm/mref_upper_err,E_cm/mref_lower_err,E_cm/"
                "mref_sym_err,"
             << "(q/mref)^2,(q/mref)^2_upper_err,(q/mref)^2_lower_err,(q/"
                "mref)^2_sym_err,"
             << "q/mref_cot_delta,q/mref_cot_delta_upper_err,q/"
                "mref_cot_delta_lower_err,"
             << "q/mref_cot_delta_sym_err" << endl;
      }

      // Get energies for this specific block using our proper mapping
      const vector<uint>& energy_indices = block_energy_indices[blocknum];

      // Process each energy in this block (unified logic for both full and
      // resampled modes)
      for (uint idx : energy_indices) {
        const RVector& energy_values =
            m_obs->getFullAndSamplingValues(labenergies[idx]);

        // Convert energy to ratio if needed (following chisq_detres.cc pattern)
        const RVector* labenergy = &energy_values;
        if (!energy_ratios) {
          KBObsInfo scalekey(mcens, MCObsInfo("KBScale"));
          const RVector& atrefmass = m_obs->getFullAndSamplingValues(scalekey);
          RVector labenergyratio(energy_values);
          if (labenergyratio.size() != atrefmass.size())
            throw(std::runtime_error(
                "Size mismatch while forming energy ratios"));
          for (uint kk = 0; kk < labenergyratio.size(); ++kk)
            labenergyratio[kk] /= atrefmass[kk];
          labenergy = &(m_obs->putFullAndSamplings(labenergies[idx],
                                                   labenergyratio, true));
        }

        // Storage for resampling results
        vector<double> elab_values, ecm_values, qsqr_values, qcot_values;
        uint nsamples = (outmode == "full") ? 1 : nsamp + 1;

        // Process full sample and resamplings
        for (uint b = 0; b < nsamples; ++b) {
          const double& Elab_ref = (*labenergy)[b]; // Elab/mref

          // Set masses for this sample
          bqptr->setRefMassL(mrefL[b]);
          bqptr->setMassesOverRef(0, (*particlemass1)[b], (*particlemass2)[b]);

          // Convert to center-of-mass energy
          double Ecm_ref = bqptr->getEcmOverMrefFromElab(Elab_ref);

          // Get center-of-mass momentum for single channel
          RVector qcmsq_over_mrefsq;
          bqptr->getQcmsqOverMrefsqFromElab(Elab_ref, qcmsq_over_mrefsq);
          double qsqr_over_mrefsq = qcmsq_over_mrefsq[0];

          // For our special case, q_cm/mref * cot(delta_0) = B_{00}
          ComplexHermitianMatrix B;
          bqptr->getBoxMatrixFromEcm(Ecm_ref, B);
          double q_cot_delta = B.get(0, 0).real(); // Take real part

          // Store results
          elab_values.push_back(Elab_ref);
          ecm_values.push_back(Ecm_ref);
          qsqr_values.push_back(qsqr_over_mrefsq);
          qcot_values.push_back(q_cot_delta);
        }

        // Output results based on mode
        if (outmode == "full") {
          fout << elab_values[0] << "," << ecm_values[0] << ","
               << qsqr_values[0] << "," << qcot_values[0] << endl;
        } else { // resampled mode
          // Calculate error estimates using jackknife or bootstrap
          MCEstimate elab_est, ecm_est, qsqr_est, qcot_est;

          // Convert vectors to RVector for analysis
          RVector elab_rvec(elab_values), ecm_rvec(ecm_values);
          RVector qsqr_rvec(qsqr_values), qcot_rvec(qcot_values);

          if (m_obs->isJackknifeMode()) {
            m_obs->jack_analyze(elab_rvec, elab_est);
            m_obs->jack_analyze(ecm_rvec, ecm_est);
            m_obs->jack_analyze(qsqr_rvec, qsqr_est);
            m_obs->jack_analyze(qcot_rvec, qcot_est);
          } else {
            m_obs->boot_analyze(elab_rvec, elab_est);
            m_obs->boot_analyze(ecm_rvec, ecm_est);
            m_obs->boot_analyze(qsqr_rvec, qsqr_est);
            m_obs->boot_analyze(qcot_rvec, qcot_est);
          }

          // Output format: full_value upper_error down_error
          double elab_upper =
              elab_est.getUpperConfLimit() - elab_est.getAverageEstimate();
          double elab_lower =
              elab_est.getAverageEstimate() - elab_est.getLowerConfLimit();
          double elab_sym_err = elab_est.getSymmetricError();
          double ecm_upper =
              ecm_est.getUpperConfLimit() - ecm_est.getAverageEstimate();
          double ecm_lower =
              ecm_est.getAverageEstimate() - ecm_est.getLowerConfLimit();
          double ecm_sym_err = ecm_est.getSymmetricError();
          double qsqr_upper =
              qsqr_est.getUpperConfLimit() - qsqr_est.getAverageEstimate();
          double qsqr_lower =
              qsqr_est.getAverageEstimate() - qsqr_est.getLowerConfLimit();
          double qsqr_sym_err = qsqr_est.getSymmetricError();
          double qcot_upper =
              qcot_est.getUpperConfLimit() - qcot_est.getAverageEstimate();
          double qcot_lower =
              qcot_est.getAverageEstimate() - qcot_est.getLowerConfLimit();
          double qcot_sym_err = qcot_est.getSymmetricError();

          fout << elab_values[0] << "," << elab_upper << "," << elab_lower
               << "," << elab_sym_err << "," << ecm_values[0] << ","
               << ecm_upper << "," << ecm_lower << "," << ecm_sym_err << ","
               << qsqr_values[0] << "," << qsqr_upper << "," << qsqr_lower
               << "," << qsqr_sym_err << "," << qcot_values[0] << ","
               << qcot_upper << "," << qcot_lower << "," << qcot_sym_err
               << endl;
        }
      }

      fout.close();
    }
    // Clean up BoxQuantization objects
    for (uint k = 0; k < BQ.size(); ++k)
      delete BQ[k];
    m_obs->clearSamplings();

  } catch (const std::exception& xp) {
    string msg("doSingleChannel failed: ");
    msg += xp.what();
    cout << msg << endl;
    throw(std::runtime_error(msg));
  }
}
