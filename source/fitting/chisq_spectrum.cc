/**
 * @file chisq_spectrum.cc
 * @brief Implementation of spectrum fitting methodology for K-matrix parameter
 * estimation
 *
 * This file contains the implementation of the SpectrumFit class, which
 * provides a comprehensive framework for fitting scattering parameters (via the
 * K-matrix) by directly comparing to the measured energy spectrum using the
 * Lüscher method.
 *
 * The spectrum fitting approach differs from determinant residual fitting by:
 * 1. Solving the Omega function Ω(Ecm) = 0 to predict energy shifts from
 * non-interacting levels
 * 2. Computing χ² residuals between predicted and observed energy shifts
 * 3. Including prior contributions for non-dElab observables
 *
 * Performance optimizations implemented:
 * - Pre-allocated temporary vectors to avoid repeated memory allocations
 * - Cached frequently accessed values to reduce computational overhead
 * - Cache-friendly memory access patterns in critical loops
 * - Minimized function call overhead in hot paths
 *
 * @author KBfit Development Team
 * @date 2024
 */

#include "chisq_spectrum.h"
#include "task_utils.h"
#include <limits>

// **************************************************************************
// *                         SPECTRUM FITTING METHODOLOGY                   *
// **************************************************************************
// *                                                                        *
// *    The class "SpectrumFit", derived from the base class "ChiSquare",   *
// *    provides a comprehensive framework for fitting scattering           *
// *    parameters (via the K-matrix) by directly comparing to the          *
// *    measured energy spectrum. It evaluates the chi^2 value by           *
// *    simultaneously treating both the K-matrix parameters and the        *
// *    underlying lattice QCD parameters (particle masses, lattice scale,  *
// *    anisotropy) as fit parameters.                                      *
// *                                                                        *
// *    The total chi-squared is constructed from two types of residuals,   *
// *    fully accounting for their correlations:                            *
// *                                                                        *
// *    1. Energy Residuals: (E_lab_shift_measured - E_lab_shift_predicted) *
// *       The predicted energy shifts are found by numerically finding     *
// *       the roots of the quantization condition (e.g., det(Kinv-B)=0)    *
// *       to predict energy shifts from non-interacting energies.         *
// *                                                                        *
// *    2. Prior Residuals: (param_prior - param_fit_value)                 *
// *       For each lattice QCD parameter that is not held fixed, this term *
// *       constrains the fit parameter to its prior value determined from  *
// *       Monte Carlo measurements.                                        *
// *                                                                        *
// *    The full covariance matrix between all measured quantities          *
// *    (energy shifts and lattice parameters) is used.                     *
// *                                                                        *
// *    XML format for chi-square fitting (Spectrum Method):                *
// *                                                                        *
// *    <SpectrumFit>                                                       *
// *                                                                        *
// *      <OmegaMu>8.0</OmegaMu>  (optional, for Stilde-based QCs)           *
// *                                                                        *
// *      <Verbose/>  (optional)                                            *
// *                                                                        *
// *      <KtildeMatrix> or <KtildeMatrixInverse>...</KtildeMatrix>         *
// *                                                                        *
// *      <DefaultEnergyFormat>...</DefaultEnergyFormat>                    *
// *                                                                        *
// *      <RootFinder> (Mandatory)                                          *
// *        Configures the adaptive bracketing root-finding algorithm.      *
// *         <MaxIterations>100</MaxIterations> (optional, default 100)     *
// *         <InitialBracketFactor>1.2</InitialBracketFactor> (opt, def 1.2)*
// *              Factor to expand search bracket if root is not found.     *
// *         <Tolerance>1e-9</Tolerance> (optional, default 1e-9)            *
// *              The tolerance for the root value itself (the y-value).    *
// *         <XTolerance>1e-9</XTolerance> (optional, default 1e-9)          *
// *              The tolerance for the energy (the x-value).               *
// *      </RootFinder>                                                     *
// *                                                                        *
// *      <MCEnsembleParameters>...</MCEnsembleParameters>                  *
// *                                                                        *
// *      <KBBlock>...</KBBlock>                                            *
// *        ...                                                             *
// *        <LabFrameEnergyShift>                                           *
// *          <MCObs>...</MCObs>                                            *
// *          <NonInteractingPair>pi(1)pi(0)</NonInteractingPair>         *
// *        </LabFrameEnergyShift>                                          *
// *        <LabFrameEnergyMin>_val_</LabFrameEnergyMin> (Mandatory)         *
// *        <LabFrameEnergyMax>_val_</LabFrameEnergyMax> (Mandatory)         *
// *                         OR                                               *
// *        <CMFrameEnergyAutoRangeMargin>_val_</CMFrameEnergyAutoRangeMargin>
// (optional) *
// *        ...                                                             *
// *                                                                        *
// *      <KBObservables> ... </KBObservables>                               *
// *                                                                        *
// *    </SpectrumFit>                                                      *
// *                                                                        *
// **************************************************************************

using namespace std;

SpectrumFit::SpectrumFit()
    : ChiSquare(), KBOH(nullptr), Kmat(nullptr), Kinv(nullptr),
      n_kmat_params(0), n_decay_channels(0), omega_mu(-1.0) {
  // Initialize to default/empty state - will be populated by clone method
}

SpectrumFit::SpectrumFit(XMLHandler& xmlin, KBObsHandler* kboh,
                         XMLHandler& xmlout, const string& outfile_stub) {
  KBOH = kboh;
  Kmat = 0;
  Kinv = 0;
  try {
    XMLHandler xmlf(xmlin, "SpectrumFit");

    bool verbose = false;
    if (xml_tag_count(xmlf, "Verbose") > 0)
      verbose = true;
    stringstream logger;
    logger.precision(12);

    //  get sampling mode and names of input files
    XMLHandler xmlr(xmlin, "KBObservables");
    MCSamplingInfo sampinfo(xmlr);
    if (sampinfo != KBOH->getSamplingInfo()) {
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

    //  first set up the K or K^(-1) matrix // TODO: allow use of K or Kinv for
    //  cayley transfored QCs
    int k1 = xmlf.count_among_children("KtildeMatrix");
    int k2 = xmlf.count_among_children("KtildeMatrixInverse");
    if ((k1 + k2) != 1)
      throw(std::invalid_argument(
          "A single KtildeMatrix or KtildeMatrixInverse tag must be present"));
    const vector<DecayChannelInfo>* dcptr;
    if (k1 == 1) {
      Kmat = new KtildeMatrixCalculator(xmlf, true);
      n_kmat_params = Kmat->getNumberOfParameters();
      n_decay_channels = Kmat->getNumberOfDecayChannels();
      dcptr = &(Kmat->getDecayChannelInfos());
    } else {
      Kinv = new KtildeInverseCalculator(xmlf, true);
      n_kmat_params = Kinv->getNumberOfParameters();
      n_decay_channels = Kinv->getNumberOfDecayChannels();
      dcptr = &(Kinv->getDecayChannelInfos());
    }

    // get set of particle names
    set<string> pnames; // same for each ensemble
    are_decay_channels_identical.resize(n_decay_channels);
    for (uint dc = 0; dc < n_decay_channels; ++dc) {
      const DecayChannelInfo& dci = (*dcptr)[dc];
      pnames.insert(dci.getParticle1Name());
      pnames.insert(dci.getParticle2Name());
      are_decay_channels_identical[dc] =
          (dci.getParticle1Name() == dci.getParticle2Name());
    }

    //  assign mu  (negative value means use determinant itself)
    omega_mu = -1.0;
    xmlreadif(xmlin, "OmegaMu", omega_mu, "SpectrumFit");
    logger << "Omega mu = " << omega_mu << endl;

    //  get format for energies/mass (in terms of product with
    //  lattice time spacing or a ratio with reference mass)
    bool energy_ratios = true;
    {
      string reply;
      xmlreadif(xmlin, "DefaultEnergyFormat", reply, "SpectrumFit");
      if ((reply == "reference_ratio") || (reply.empty()))
        energy_ratios = true;
      else if (reply == "time_spacing_product")
        energy_ratios = false;
      else
        throw(std::invalid_argument("Invalid <DefaultEnergyFormat> tag"));
    }

    // Get root finding info and initialize energy bounds
    int root_counts = xmlin.count_among_children("RootFinder");
    if (root_counts != 1)
      throw(std::invalid_argument("There must be one only RootFinder tag"));

    XMLHandler xmlroot(xmlin, "RootFinder");
    root_finder_config = AdaptiveBracketRootFinder::makeConfigFromXML(xmlroot);

    // otherwise, default config
    logger << "RootFinder configuration: " << endl;
    logger << root_finder_config.toString() << endl;

    //  Loop over the different MC ensembles to get information
    //  about all needed parameters.  Particle masses, anisotropy,
    //  parameters should be associated with the ensemble.

    map<MCEnsembleInfo, set<MCObsInfo>> needed_keys;
    set<MCEnsembleInfo> ensembles;
    map<MCEnsembleInfo, MCObsInfo> ref_at_mass;
    map<MCEnsembleInfo, MCObsInfo> anisotropy;
    map<MCEnsembleInfo, map<string, MCObsInfo>> particle_masses;
    map<KBObsInfo, double> fixed_values;

    list<XMLHandler> xmlen = xmlf.find("MCEnsembleParameters");
    for (list<XMLHandler>::iterator it = xmlen.begin(); it != xmlen.end();
         ++it) {
      MCEnsembleInfo mcens(*it);
      EnsembleFitData this_ensemble_data(mcens);

      if (ensembles.find(mcens) != ensembles.end())
        throw(std::invalid_argument("Duplicate MC ensemble"));
      ensembles.insert(mcens);
      set<MCObsInfo>& kset = needed_keys[mcens];
      string pname;
      MCObsInfo rkey;
      // get lattice spacing times reference scale info
      read_obs(*it, "ReferenceMassTimeSpacingProduct", false, rkey, kset, pname,
               mcens, fixed_values);
      ref_at_mass.insert(make_pair(mcens, rkey));
      if (!energy_ratios) {
        kset.erase(rkey);
      } // remove ref_at_mass from kset
        // get anisotropy info
      if (xml_tag_count(*it, "LatticeAnisotropy") == 1) {
        read_obs(*it, "LatticeAnisotropy", false, rkey, kset, pname, mcens,
                 fixed_values);
        anisotropy.insert(make_pair(mcens, rkey));
      } else if (xml_tag_count(*it, "LatticeAnisotropy") > 1) {
        throw(std::invalid_argument(
            "Multiple LatticeAnisotropy tags cannot be present"));
      }
      // get particle mass infos
      map<string, MCObsInfo> pmap;
      list<XMLHandler> xmlp = it->find("ParticleMass");
      for (list<XMLHandler>::iterator pt = xmlp.begin(); pt != xmlp.end();
           ++pt) {
        read_obs(*pt, "ParticleMass", true, rkey, kset, pname, mcens,
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

    //  connect files for input

    KBOH->connectSamplingFiles(sampfiles, needed_keys, verbose);
    logger << KBOH->getCurrentLog().str();
    KBOH->clearLog();

    //  Loop over the KB quantization blocks to get the shift
    //  energies infos and samplings

    vector<MCEnsembleInfo> blockens;
    map<MCEnsembleInfo, uint> ensemble_idmap;

    list<XMLHandler> xmlkb = xmlf.find("KBBlock");
    uint blockcount = 0;
    uint ensemblecount = 0;
    for (list<XMLHandler>::iterator it = xmlkb.begin(); it != xmlkb.end();
         ++it) {
      BoxQuantization* bqptr = new BoxQuantization(*it, Kmat, Kinv);
      MCEnsembleInfo mcens(*it);
      if (ensembles.find(*it) == ensembles.end())
        throw(
            std::invalid_argument("KBBlock associated with unknown ensemble"));
      blockens.push_back(mcens);
      if (ensemble_idmap.find(mcens) == ensemble_idmap.end()) {
        ensemble_idmap.insert(make_pair(mcens, ensemblecount));
        ensemble_fit_data.resize(ensemblecount + 1);
        ensemble_fit_data[ensemblecount].ensemble_info = mcens;
        ensemble_fit_data[ensemblecount].n_blocks =
            0; // Will be incremented below
        ensemblecount++;
        // assumes that the xml would be ordered by ensemble, then blocks
      }
      set<MCObsInfo>& kset = needed_keys[mcens];
      uint nres = 0;
      MCObsInfo rkey;
      list<XMLHandler> xmlee = it->find("LabFrameEnergyShift");
      for (auto& eet : xmlee) {
        read_obs(eet, "LabFrameEnergyShift", rkey, kset);

        // Read the NonInteractingPair information
        string ni_pair_str;
        try {
          xmlread(eet, "NonInteractingPair", ni_pair_str,
                  "LabFrameEnergyShift");
        } catch (const std::exception& xp) {
          throw(std::invalid_argument("NonInteractingPair tag must be present "
                                      "in LabFrameEnergyShift: " +
                                      string(xp.what())));
        }
        NonInteractingPair NI_pair =
            parseNonInteractingPair(ni_pair_str, *dcptr);
        ensemble_fit_data[ensemble_idmap[mcens]]
            .non_interacting_pairs.push_back(NI_pair);

        KBObsInfo this_energy_key(mcens, rkey);
        ensemble_fit_data[ensemble_idmap[mcens]].dElab_samples.push_back(
            KBOH->getFullAndSamplingValues(this_energy_key));
        ensemble_fit_data[ensemble_idmap[mcens]].energy_obs_infos.push_back(
            rkey);
        nres++;
      }
      if (nres == 0)
        throw(std::invalid_argument(
            "No energies available in at least one block"));

      // set the energy bounds for this block
      double Ecm_min = 0.0, Ecm_max = 0.0;
      bool cm_frame_min_given =
          xmlreadifchild(*it, "CMFrameEnergyMin", Ecm_min);
      bool cm_frame_max_given =
          xmlreadifchild(*it, "CMFrameEnergyMax", Ecm_max);

      if (!(cm_frame_min_given && cm_frame_max_given)) {
        double auto_margin = 0.0;
        bool auto_bounds =
            xmlreadifchild(*it, "AutoEcmBoundsMargin", auto_margin);
        if (!auto_bounds)
          throw(std::invalid_argument("CMFrameEnergyMin/Max missing and "
                                      "AutoEcmBoundsMargin not provided"));

        size_t start_idx = ensemble_fit_data[ensemble_idmap[mcens]]
                               .non_interacting_pairs.size() -
                           nres;
        const auto& pairs =
            ensemble_fit_data[ensemble_idmap[mcens]].non_interacting_pairs;
        double emin = std::numeric_limits<double>::infinity();
        double emax = 0.0;
        for (uint jj = 0; jj < nres; ++jj) {
          const NonInteractingPair& ni = pairs[start_idx + jj];
          const EcmTransform& et =
              bqptr->getDecayChannelEcmTransform(ni.decay_channel_idx);
          double efree = et.getFreeTwoParticleEnergyInEcm(ni.d1_sqr, ni.d2_sqr);
          emin = std::min(emin, efree);
          emax = std::max(emax, efree);
        }
        Ecm_min = emin - auto_margin;
        Ecm_max = emax + auto_margin;
      }
      ensemble_fit_data[ensemble_idmap[mcens]]
          .Ecm_bounds_per_block.emplace_back(Ecm_min, Ecm_max);
      ensemble_fit_data[ensemble_idmap[mcens]].n_energies_per_block.push_back(
          nres);
      ensemble_fit_data[ensemble_idmap[mcens]].BQ_blocks.push_back(bqptr);
      ensemble_fit_data[ensemble_idmap[mcens]].n_blocks++;
      blockcount++;
    }
    if (blockcount == 0) {
      throw(std::runtime_error("No data to analyze"));
    }

    // Calculate total number of fit parameters and residuals
    uint total_fit_params = n_kmat_params; // K-matrix parameters
    uint total_residuals = 0;
    uint numenergies = 0;

    for (uint k = 0; k < ensemble_fit_data.size(); ++k) {
      // mref is always a parameter
      total_fit_params++;
      total_residuals++;

      // anisotropy if not fixed
      if (!ensemble_fit_data[k].is_anisotropy_fixed) {
        total_fit_params++;
        total_residuals++;
      }

      uint non_fixed_masses =
          ensemble_fit_data[k].mass_samples.size(); // Only non-fixed masses
      total_fit_params += non_fixed_masses;
      total_residuals += non_fixed_masses;

      // Count energy residuals (these are observations, don't contribute to fit
      // parameters)
      for (uint n_energies : ensemble_fit_data[k].n_energies_per_block) {
        numenergies += n_energies;
        total_residuals += n_energies;
      }
    }

    // get QuantizationCondition
    string qctype;
    xmlreadif(xmlf, "QuantizationCondition", qctype, "SpectrumFit");
    if (qctype.empty()) {
      throw(std::invalid_argument("QuantizationCondition tag must be present"));
    }
    BoxQuantization* bqptr_dummy = ensemble_fit_data[0].BQ_blocks[0];
    try {
      qctype_enum = bqptr_dummy->getQuantCondTypeFromString(qctype).value();
      if (qctype_enum == BoxQuantization::KtildeB) {
        if (k2 == 1) {
          throw(std::invalid_argument(
              "KtildeMatrixInverse cannot be used with KtildeB"));
        }
      }
      if (qctype_enum == BoxQuantization::KtildeinvB) {
        if (k1 == 1) {
          throw(std::invalid_argument(
              "KtildeMatrix cannot be used with KtildeinvB"));
        }
      }
    } catch (std::bad_optional_access&) {
      throw(std::invalid_argument("Invalid QuantizationCondition tag"));
    }

    //  get number of resamplings.  For bootstrap, all ensembles
    //  must have same bootstrap parameters.  For jackknife,
    //  all ensembles must have the SAME number of bins.

    uint nsamplings = KBOH->getNumberOfResamplings();
    if (KBOH->isJackknifeMode()) {
      for (set<MCEnsembleInfo>::const_iterator et = ensembles.begin();
           et != ensembles.end(); ++et) {
        if (KBOH->getNumberOfBins(*et) != nsamplings) {
          throw(std::invalid_argument(
              "All ensembles must have same number of jackknife bins"));
        }
      }
    }

    if (nsamplings == 0)
      throw(std::runtime_error(
          "Samplings numbers do not match for all ensembles"));

    //  insert all of the fixed values into KBOH

    for (map<KBObsInfo, double>::const_iterator fx = fixed_values.begin();
         fx != fixed_values.end(); ++fx) {
      KBOH->putFixedValue(fx->first, fx->second, KBOH->getNumberOfResamplings(),
                          true);
    }

    //  read the reference mass time-spacing products for each ensemble into
    //  memory; evaluate reference lengths and particle masses for each ensemble
    for (const auto& et : ref_at_mass) {
      const MCEnsembleInfo& mcens = et.first;
      // Ensure that the ensemble from ref_at_mass is actually used in a
      // KBBlock, otherwise, we have no need for it.
      if (ensemble_idmap.find(mcens) == ensemble_idmap.end()) {
        continue;
      }
      uint current_ensemble_id = ensemble_idmap.at(mcens);
      uint nsamp = KBOH->getNumberOfResamplings();
      KBObsInfo atrefmasskey(mcens, et.second);
      const RVector& atrefmass0 = KBOH->getFullAndSamplingValues(atrefmasskey);
      if (atrefmass0.size() != (nsamp + 1))
        throw(std::runtime_error("Resampling size mismatch in KBfit"));
      KBObsInfo scalekey(mcens, MCObsInfo("KBScale")); // ref mass
      KBOH->putFullAndSamplings(scalekey, atrefmass0, true);
      KBOH->eraseSamplings(atrefmasskey);

      // get the reference length

      uint Llat = mcens.getLatticeXExtent();
      if ((mcens.getLatticeYExtent() != Llat) ||
          (mcens.getLatticeZExtent() != Llat))
        throw(
            std::runtime_error("KBfit only works for LxLxL spatial lattices"));
      const RVector& atrefmass = KBOH->getFullAndSamplingValues(scalekey);
      if (atrefmass.size() != (nsamp + 1))
        throw(std::runtime_error("Resampling size mismatch in KBfit"));

      // Store lattice extent for this ensemble
      ensemble_fit_data[current_ensemble_id].Llat = Llat;

      // Convert energy values to reference mass ratios if in
      // time_spacing_product mode
      if (!energy_ratios) {
        EnsembleFitData& ens_data = ensemble_fit_data[current_ensemble_id];
        for (uint energy_idx = 0; energy_idx < ens_data.dElab_samples.size();
             ++energy_idx) {
          RVector& energy_samples = ens_data.dElab_samples[energy_idx];
          if (energy_samples.size() != atrefmass.size())
            throw(std::runtime_error(
                "Size mismatch while forming energy ratios"));
          for (uint kk = 0; kk < energy_samples.size(); ++kk) {
            energy_samples[kk] /= atrefmass[kk];
          }
        }
      }

      // Always store mref as an observable
      ensemble_fit_data[current_ensemble_id].mref_samples.resize(nsamp + 1);
      ensemble_fit_data[current_ensemble_id].mref_samples = atrefmass;
      ensemble_fit_data[current_ensemble_id].mref_prior =
          MCObsInfo("KBScale", current_ensemble_id);

      // Handle anisotropy
      if (anisotropy.find(mcens) == anisotropy.end()) {
        // Isotropic case - anisotropy is fixed to 1.0
        ensemble_fit_data[current_ensemble_id].is_anisotropy_fixed = true;
        ensemble_fit_data[current_ensemble_id].fixed_anisotropy_value = 1.0;
      } else {
        // Anisotropic case - check if it's fixed
        KBObsInfo anisotropykey(mcens, anisotropy[mcens]);
        bool anisotropy_fixed =
            (fixed_values.find(anisotropykey) != fixed_values.end());
        ensemble_fit_data[current_ensemble_id].is_anisotropy_fixed =
            anisotropy_fixed;

        if (anisotropy_fixed) {
          ensemble_fit_data[current_ensemble_id].fixed_anisotropy_value =
              fixed_values[anisotropykey];
        } else {
          const RVector& anisotropy_samples =
              KBOH->getFullAndSamplingValues(anisotropykey);
          if (anisotropy_samples.size() != (nsamp + 1))
            throw(std::runtime_error("Resampling size mismatch in KBfit"));
          ensemble_fit_data[current_ensemble_id].anisotropy_samples =
              anisotropy_samples;
          ensemble_fit_data[current_ensemble_id].anisotropy_prior =
              MCObsInfo("LatticeAnisotropy", current_ensemble_id);
        }
      }

      // Note: Length is always calculated from mref * Llat * anisotropy

      //  get particle masses for each decay channel (form ratios if
      //  energy_ratio false) Only add non-fixed masses to observables, store
      //  fixed values separately

      map<string, MCObsInfo>& pmap = particle_masses[mcens];

      // Initialize fixed mass values vector (2 particles per decay channel max)
      ensemble_fit_data[current_ensemble_id].fixed_mass_values.resize(
          n_decay_channels * 2, 0.0);
      ensemble_fit_data[current_ensemble_id].is_mass_fixed.resize(
          n_decay_channels * 2, false);

      // Process masses for each decay channel
      for (uint dc = 0; dc < n_decay_channels; ++dc) {
        const DecayChannelInfo& dci = (*dcptr)[dc];

        // First particle in decay channel
        string particle1_name = dci.getParticle1Name();
        MCObsInfo& particle1_key = pmap[particle1_name];
        particle1_key.resetObsIndex(dc * 2);
        KBObsInfo mass1key(mcens, particle1_key);

        bool mass1_fixed = (fixed_values.find(mass1key) != fixed_values.end());
        uint mass1_idx = dc * 2;

        ensemble_fit_data[current_ensemble_id].is_mass_fixed[mass1_idx] =
            mass1_fixed;

        if (mass1_fixed) {
          ensemble_fit_data[current_ensemble_id].fixed_mass_values[mass1_idx] =
              fixed_values[mass1key];
        } else {
          // Add to observables
          const RVector& mass1 = KBOH->getFullAndSamplingValues(mass1key);
          ensemble_fit_data[current_ensemble_id].mass_priors.push_back(
              particle1_key);
          if ((!energy_ratios)) {
            RVector massratio(mass1);
            if (massratio.size() != atrefmass.size())
              throw(std::runtime_error(
                  "Size mismatch while forming mass ratios"));
            for (uint kk = 0; kk < massratio.size(); ++kk)
              massratio[kk] /= atrefmass[kk];
            KBOH->putFullAndSamplings(mass1key, massratio, true);
            ensemble_fit_data[current_ensemble_id].mass_samples.push_back(
                massratio);
          } else {
            ensemble_fit_data[current_ensemble_id].mass_samples.push_back(
                mass1);
          }
        }

        // Second particle in decay channel (only if different from first)
        if (!are_decay_channels_identical[dc]) {
          string particle2_name = dci.getParticle2Name();
          MCObsInfo& particle2_key = pmap[particle2_name];
          particle2_key.resetObsIndex(dc * 2 + 1);
          KBObsInfo mass2key(mcens, particle2_key);

          bool mass2_fixed =
              (fixed_values.find(mass2key) != fixed_values.end());
          uint mass2_idx = dc * 2 + 1;

          ensemble_fit_data[current_ensemble_id].is_mass_fixed[mass2_idx] =
              mass2_fixed;

          if (mass2_fixed) {
            // Store fixed value
            ensemble_fit_data[current_ensemble_id]
                .fixed_mass_values[mass2_idx] = fixed_values[mass2key];
          } else {
            // Add to observables
            const RVector& mass2 = KBOH->getFullAndSamplingValues(mass2key);
            ensemble_fit_data[current_ensemble_id].mass_priors.push_back(
                particle2_key);
            if ((!energy_ratios)) {
              RVector massratio(mass2);
              if (massratio.size() != atrefmass.size())
                throw(std::runtime_error(
                    "Size mismatch while forming mass ratios"));
              for (uint kk = 0; kk < massratio.size(); ++kk)
                massratio[kk] /= atrefmass[kk];
              KBOH->putFullAndSamplings(mass2key, massratio, true);
              ensemble_fit_data[current_ensemble_id].mass_samples.push_back(
                  massratio);
            } else {
              ensemble_fit_data[current_ensemble_id].mass_samples.push_back(
                  mass2);
            }
          }
        }
      }
    }
    initialize_base(total_fit_params, total_residuals, nsamplings);

    // Initialize decay channel masses vector
    decay_channel_masses.resize(n_decay_channels);

    // Pre-allocate temporary vectors for performance (estimate max size)
    uint max_energies_per_block = 0;
    for (const auto& ens_data : ensemble_fit_data) {
      for (uint n_energies : ens_data.n_energies_per_block) {
        max_energies_per_block = std::max(max_energies_per_block, n_energies);
      }
    }
    energy_shift_predictions.reserve(max_energies_per_block);
    shift_obs_w_NIs.reserve(max_energies_per_block);
    fn_calls.reserve(max_energies_per_block);

    // The SpectrumFit is designed to use a fixed covariance matrix calculated
    // once.
    this->initializeInvCovCholesky();

    KBOH->clearSamplings(); // Clear sampling data from memory after initial use
    xmlformat("SpectrumFitConstruction", logger.str(), xmlout);

  } catch (const std::exception& xp) {
    string msg("Construction of SpectrumFit failed: ");
    msg += xp.what();
    cout << msg << endl;
    clear();
    throw(std::runtime_error(msg));
  }
}

SpectrumFit::~SpectrumFit() { clear(); }

void SpectrumFit::clear() {
  for (uint ens = 0; ens < ensemble_fit_data.size(); ++ens) {
    for (uint k = 0; k < ensemble_fit_data[ens].BQ_blocks.size(); ++k)
      delete ensemble_fit_data[ens].BQ_blocks[k];
    ensemble_fit_data[ens].BQ_blocks.clear();
  }
  delete Kmat;
  delete Kinv;
  ensemble_fit_data.clear();

  // Clear reusable vectors
  energy_shift_predictions.clear();
  fn_calls.clear();
  decay_channel_masses.clear();
}

std::unique_ptr<SpectrumFit> SpectrumFit::clone(KBObsHandler* new_kboh) const {
  // Create a new instance using the private default constructor
  auto cloned = std::unique_ptr<SpectrumFit>(new SpectrumFit());

  // Copy base class members
  cloned->nresiduals = this->nresiduals;
  cloned->nfitparams = this->nfitparams;
  cloned->inv_cov_cholesky = this->inv_cov_cholesky;
  cloned->residuals = this->residuals;
  cloned->nresamplings = this->nresamplings;
  cloned->resampling_index = this->resampling_index;
  cloned->qctype_enum = this->qctype_enum;

  // Copy or assign KBObsHandler pointer (shallow copy - external object)
  cloned->KBOH = (new_kboh != nullptr) ? new_kboh : this->KBOH;

  // Copy configuration data (value types - deep copied automatically)
  cloned->root_finder_config = this->root_finder_config;
  cloned->n_kmat_params = this->n_kmat_params;
  cloned->n_decay_channels = this->n_decay_channels;
  cloned->omega_mu = this->omega_mu;
  cloned->are_decay_channels_identical = this->are_decay_channels_identical;

  // Deep copy K-matrix calculators first (needed for BoxQuantization cloning)
  if (this->Kmat != nullptr) {
    // Use the new clone method instead of copy constructor
    cloned->Kmat = this->Kmat->clone().release();
    cloned->Kinv = nullptr;
  } else if (this->Kinv != nullptr) {
    // Use the new clone method instead of copy constructor
    cloned->Kinv = this->Kinv->clone().release();
    cloned->Kmat = nullptr;
  } else {
    cloned->Kmat = nullptr;
    cloned->Kinv = nullptr;
  }

  // Deep copy ensemble fit data
  cloned->ensemble_fit_data.clear();
  cloned->ensemble_fit_data.reserve(this->ensemble_fit_data.size());

  for (const auto& ens_data : this->ensemble_fit_data) {
    EnsembleFitData cloned_ens_data =
        ens_data; // Copy most members via default copy semantics

    // Deep copy BoxQuantization pointers within each ensemble
    cloned_ens_data.BQ_blocks.clear();
    cloned_ens_data.BQ_blocks.reserve(ens_data.BQ_blocks.size());
    for (const auto* bq : ens_data.BQ_blocks) {
      if (bq != nullptr) {
        // Use the new clone method with the cloned K-matrix calculators
        std::unique_ptr<BoxQuantization> cloned_bq =
            bq->clone(cloned->Kmat, cloned->Kinv);
        cloned_ens_data.BQ_blocks.push_back(cloned_bq.release());
      } else {
        cloned_ens_data.BQ_blocks.push_back(nullptr);
      }
    }

    cloned->ensemble_fit_data.push_back(std::move(cloned_ens_data));
  }

  // Initialize temporary vectors (these are mutable working space - start
  // empty)
  cloned->energy_shift_predictions.clear();
  cloned->fn_calls.clear();
  cloned->decay_channel_masses.resize(this->decay_channel_masses.size());

  return cloned;
}

// Need to add non-Kmatrix fit parameters below here
// (L, mass, etc.), just use mean values.
void SpectrumFit::guessInitialFitParamValues(vector<double>& fitparams,
                                             bool only_update_priors) const {
  if (!only_update_priors) {
    // 1. K-matrix parameters
    const std::vector<double>& k_fit_params =
        (Kmat ? Kmat->getParameterValues() : Kinv->getParameterValues());
    // Copy K-matrix params in one go
    std::copy(k_fit_params.begin(), k_fit_params.end(), fitparams.begin());
  }
  std::size_t offset = n_kmat_params;

  // 2. Ensemble-specific prior parameters: only non-fixed ones
  for (std::size_t e = 0; e < ensemble_fit_data.size(); ++e) {
    const EnsembleFitData& ens_data = ensemble_fit_data[e]; // Cache reference

    // mref (always present)
    fitparams[offset++] = ens_data.mref_samples[resampling_index];

    // Anisotropy (only if not fixed)
    if (!ens_data.is_anisotropy_fixed) {
      fitparams[offset++] = ens_data.anisotropy_samples[resampling_index];
    }

    // Masses (only non-fixed ones)
    for (const auto& mass_samples : ens_data.mass_samples) {
      fitparams[offset++] = mass_samples[resampling_index];
    }
  }
}

void SpectrumFit::getFitParamMCObsInfo(vector<MCObsInfo>& fitinfos) const {
  const vector<KFitParamInfo>& k_infos =
      (Kmat ? Kmat->getFitParameterInfos() : Kinv->getFitParameterInfos());

  uint offset = 0;
  // Copy K-matrix parameter infos
  for (const auto& info : k_infos) {
    fitinfos[offset++] = MCObsInfo(info.getMCObsName());
  }

  // Copy ensemble parameter infos
  for (const auto& ens_data : ensemble_fit_data) {
    // mref info (always present)
    fitinfos[offset++] = ens_data.mref_prior;

    // Anisotropy info (only if not fixed)
    if (!ens_data.is_anisotropy_fixed) {
      fitinfos[offset++] = ens_data.anisotropy_prior;
    }

    // Mass infos (only non-fixed ones)
    for (const auto& mass_prior : ens_data.mass_priors) {
      fitinfos[offset++] = mass_prior;
    }
  }
}

// might want to modify this at some point once I
// better understand its purpose
void SpectrumFit::do_output(XMLHandler& xmlout) const {
  xmlout.set_root("SpectrumFit");
  XMLHandler xmlK;
  if (Kmat != 0)
    Kmat->output(xmlK);
  else
    Kinv->output(xmlK);
  xmlout.put_child(xmlK);
  // TODO: Add any other SpectrumFit specific outputs if necessary
}

/**
 * @brief Performance-critical method for evaluating residuals and inverse
 * covariance
 * @param fitparams Current fit parameter values
 *
 * This method is the computational bottleneck of the spectrum fitting process,
 * called for every function evaluation during minimization. It performs the
 * following operations:
 *
 * 1. Parameter Extraction: Extracts K-matrix parameters, masses, and anisotropy
 * 2. Box Quantization Update: Updates BoxQuantization objects with current
 * parameters
 * 3. Root Finding: Solves Omega function Ω(Ecm) = 0 to predict energy shifts
 * from non-interacting energies
 * 4. Residual Calculation: Computes residuals between predicted and observed
 * energy shifts
 * 5. Prior Contributions: Adds prior terms for non-dElab observables
 *
 * PERFORMANCE OPTIMIZATIONS IMPLEMENTED:
 * - Pre-allocated temporary vectors to avoid repeated memory allocations
 * - Cached frequently accessed values (resampling index, length parameters)
 * - Cache-friendly memory access patterns in ensemble and block loops
 * - Minimized function call overhead through direct member access
 * - Efficient parameter passing using pointer arithmetic
 * - CENTER-OF-MASS FRAME OPTIMIZATION: Root finding now works directly in CM
 * frame to eliminate redundant Ecm ↔ Elab coordinate transformations that were
 * previously performed for every single root finding evaluation.
 *
 * The covariance matrix is fixed after initializeInvCovCholesky and not
 * recalculated.
 *
 * @note This method is called thousands of times during fitting - any
 * optimization here has significant impact on overall performance.
 */
void SpectrumFit::evalResidualsAndInvCovCholesky(
    const std::vector<double>& fitparams) {
  // OPTIMIZATION: Set K-matrix parameters directly without copying - use
  // subvector constructor
  if (Kmat) {
    Kmat->setParameterValues(std::vector<double>(
        fitparams.begin(), fitparams.begin() + n_kmat_params));
  } else {
    Kinv->setParameterValues(std::vector<double>(
        fitparams.begin(), fitparams.begin() + n_kmat_params));
  }

  // OPTIMIZATION: Use pointer arithmetic to avoid copying prior_params vector
  const double* prior_params = fitparams.data() + n_kmat_params;

  uint residual_offset = 0;    // Index into residuals array
  uint prior_param_offset = 0; // Index into prior_params array

  // OPTIMIZATION: Cache the resampling index to avoid repeated member access
  const uint current_resampling_idx = resampling_index;

  // OPTIMIZATION: Process each ensemble in cache-friendly manner
  for (uint ensemble = 0; ensemble < ensemble_fit_data.size(); ++ensemble) {
    const EnsembleFitData& ens_data =
        ensemble_fit_data[ensemble]; // OPTIMIZATION: Cache reference

    // Extract and process mref parameter (always present)
    const double mref_param = prior_params[prior_param_offset];
    residuals[residual_offset] =
        ens_data.mref_samples[current_resampling_idx] - mref_param;
    residual_offset++;
    prior_param_offset++;

    // Extract and process anisotropy parameter
    double anisotropy_param;
    if (ens_data.is_anisotropy_fixed) {
      anisotropy_param = ens_data.fixed_anisotropy_value;
    } else {
      anisotropy_param = prior_params[prior_param_offset];
      residuals[residual_offset] =
          ens_data.anisotropy_samples[current_resampling_idx] -
          anisotropy_param;
      residual_offset++;
      prior_param_offset++;
    }

    // OPTIMIZATION: Calculate length parameter once and cache it
    const double length_param =
        mref_param * double(ens_data.Llat) * anisotropy_param;

    // OPTIMIZATION: Set length parameter for all blocks in a single pass
    for (uint block = 0; block < ens_data.n_blocks; ++block) {
      ens_data.BQ_blocks[block]->setRefMassL(length_param);
    }

    uint mass_sample_idx = 0;
    for (uint decay_channel = 0; decay_channel < n_decay_channels;
         ++decay_channel) {
      const uint mass_base_idx = decay_channel * 2;
      double mass1, mass2;

      // Get mass1 (first particle in decay channel)
      if (ens_data.is_mass_fixed[mass_base_idx]) {
        mass1 = ens_data.fixed_mass_values[mass_base_idx];
      } else {
        mass1 = prior_params[prior_param_offset];
        residuals[residual_offset] =
            ens_data.mass_samples[mass_sample_idx][current_resampling_idx] -
            mass1;
        residual_offset++;
        prior_param_offset++;
        mass_sample_idx++;
      }

      // Get mass2 (second particle in decay channel)
      if (are_decay_channels_identical[decay_channel]) {
        mass2 = mass1;
      } else {
        const uint mass2_idx = mass_base_idx + 1;
        if (ens_data.is_mass_fixed[mass2_idx]) {
          mass2 = ens_data.fixed_mass_values[mass2_idx];
        } else {
          mass2 = prior_params[prior_param_offset];
          residuals[residual_offset] =
              ens_data.mass_samples[mass_sample_idx][current_resampling_idx] -
              mass2;
          residual_offset++;
          prior_param_offset++;
          mass_sample_idx++;
        }
      }

      // OPTIMIZATION: Store masses for efficient access in energy loop
      decay_channel_masses[decay_channel] = std::make_pair(mass1, mass2);

      // OPTIMIZATION: Set masses for all blocks in a single pass
      for (uint block = 0; block < ens_data.n_blocks; ++block) {
        ens_data.BQ_blocks[block]->setMassesOverRef(decay_channel, mass1,
                                                    mass2);
      }
    }

    // CRITICAL PERFORMANCE SECTION: Calculate energy residuals
    uint energy_offset = 0;
    for (uint block = 0; block < ens_data.n_blocks; ++block) {
      BoxQuantization* this_block_bq = ens_data.BQ_blocks[block];
      const uint n_energies = ens_data.n_energies_per_block[block];

      // OPTIMIZATION: Cache energy bounds for this block - work directly in CM
      // frame Energy bounds are stored in CM frame to avoid repeated
      // transformations
      const auto& bounds = ens_data.Ecm_bounds_per_block[block];
      const double Ecm_min = bounds.first;
      const double Ecm_max = bounds.second;

      // OPTIMIZATION: Reuse pre-allocated vectors (avoid repeated allocation)
      energy_shift_predictions.clear();
      energy_shift_predictions.reserve(n_energies);
      shift_obs_w_NIs.clear();
      shift_obs_w_NIs.reserve(n_energies);
      fn_calls.clear();

      // OPTIMIZATION: Cache-friendly loop over energies in this block
      // Build vector of observed energy shifts paired with non-interacting
      // pairs
      for (uint energy_index = 0; energy_index < n_energies; ++energy_index) {
        // Access non-interacting pair information for this energy level
        const uint global_energy_idx = energy_offset + energy_index;
        shift_obs_w_NIs.emplace_back(
            ens_data.dElab_samples[global_energy_idx][current_resampling_idx],
            ens_data.non_interacting_pairs[global_energy_idx]);
      }

      this_block_bq->getDeltaEnergyPredictionsOptimized(
          omega_mu, Ecm_min, Ecm_max, qctype_enum, root_finder_config,
          shift_obs_w_NIs, energy_shift_predictions, fn_calls,
          /* output_in_lab_frame = */ true);

      // Residuals are computed as: observed_shift - predicted_shift
      for (uint energy_index = 0; energy_index < n_energies; ++energy_index) {
        residuals[residual_offset++] =
            ens_data.dElab_samples[energy_offset++][current_resampling_idx] -
            energy_shift_predictions[energy_index];
      }
    }
  }
  assert(residual_offset == nresiduals);
}

/**
 * @brief Initialize inverse covariance matrix Cholesky decomposition
 *
 * This method calculates the inverse covariance matrix and its Cholesky
 * decomposition once at the beginning of the fitting process. The decomposition
 * is used by the base ChiSquare class to efficiently compute χ² values during
 * minimization.
 *
 * The covariance matrix structure is simplified by introducing priors for MC
 * observables: cov(r_i, r_j) = cov(R_i, R_j), where R is the observable and r
 * is the residual.
 *
 * The method processes observables in the following order:
 * 1. Reference mass (mref) - always present
 * 2. Anisotropy - only if not fixed
 * 3. Masses - only non-fixed ones
 * 4. Energy levels - for all ensembles and blocks
 *
 * @note This method is called once during initialization, unlike
 * evalResidualsAndInvCovCholesky which is called repeatedly during fitting.
 */
void SpectrumFit::initializeInvCovCholesky() {
  std::vector<const RVector*> obs_samples;
  std::vector<MCEnsembleInfo> obs_ensemble_infos;
  obs_samples.reserve(nresiduals);
  obs_ensemble_infos.reserve(nresiduals);

  for (uint ens = 0; ens < ensemble_fit_data.size(); ++ens) {
    const EnsembleFitData& ens_data = ensemble_fit_data[ens]; // Cache reference

    // mref (always present)
    obs_samples.push_back(&ens_data.mref_samples);
    obs_ensemble_infos.push_back(ens_data.ensemble_info);

    // anisotropy (only if not fixed)
    if (!ens_data.is_anisotropy_fixed) {
      obs_samples.push_back(&ens_data.anisotropy_samples);
      obs_ensemble_infos.push_back(ens_data.ensemble_info);
    }

    // masses (only non-fixed ones)
    uint mass_sample_idx = 0;
    for (uint dc = 0; dc < n_decay_channels; ++dc) {
      const uint mass_base_idx = dc * 2; // Cache calculation

      if (!ens_data.is_mass_fixed[mass_base_idx]) {
        obs_samples.push_back(&ens_data.mass_samples[mass_sample_idx++]);
        obs_ensemble_infos.push_back(ens_data.ensemble_info);
      }

      if (!are_decay_channels_identical[dc] &&
          !ens_data.is_mass_fixed[mass_base_idx + 1]) {
        obs_samples.push_back(&ens_data.mass_samples[mass_sample_idx++]);
        obs_ensemble_infos.push_back(ens_data.ensemble_info);
      }
    }

    // energy Delta E's - all energies from all blocks
    uint energy_offset = 0;
    for (uint block = 0; block < ens_data.n_blocks; ++block) {
      const uint nE = ens_data.n_energies_per_block[block]; // Cache value
      for (uint e = 0; e < nE; ++e) {
        obs_samples.push_back(&ens_data.dElab_samples[energy_offset++]);
        obs_ensemble_infos.push_back(ens_data.ensemble_info);
      }
    }
  }

  RealSymmetricMatrix cov(nresiduals, 0.0);
  const bool bootstrap_mode = KBOH->isBootstrapMode(); // Cache boolean

  for (uint k = 0; k < nresiduals; ++k) {
    for (uint j = 0; j <= k; ++j) {
      // Set covariance to zero if observables are from different ensembles
      if (obs_ensemble_infos[k] == obs_ensemble_infos[j]) {
        cov(k, j) =
            bootstrap_mode
                ? KBOH->boot_covariance(*(obs_samples[k]), *(obs_samples[j]))
                : KBOH->jack_covariance(*(obs_samples[k]), *(obs_samples[j]));
      }
      // else cov(k,j) remains 0.0 (already initialized)
    }
  }

  CholeskyDecomposer CD;
  CD.getCholeskyOfInverse(cov, inv_cov_cholesky);

  // DEBUG: output the covariance matrix
  // {
  //   cout << "Covariance matrix:" << endl;
  //   for (uint i = 0; i < nresiduals; ++i) {
  //     for (uint j = 0; j < nresiduals; ++j) {
  //       cout << cov(i, j) << " ";
  //     }
  //     cout << endl;
  //   }
  // }
}

NonInteractingPair SpectrumFit::parseNonInteractingPair(
    const std::string& pair_str,
    const std::vector<DecayChannelInfo>& decay_channels) const {
  NonInteractingPair result;

  // Parse format: "particle1(d1)particle2(d2)"
  // Example: "pi(1)pi(0)" or "K(2)pi(1)"
  // This represents a non-interacting pair with their individual d^2

  size_t first_paren = pair_str.find('(');
  if (first_paren == string::npos) {
    throw(std::invalid_argument(
        "Invalid NonInteractingPair format - missing first parentheses"));
  }

  string particle1_name = pair_str.substr(0, first_paren);

  size_t first_close = pair_str.find(')', first_paren);
  if (first_close == string::npos) {
    throw(std::invalid_argument("Invalid NonInteractingPair format - missing "
                                "first closing parentheses"));
  }

  string d1_str =
      pair_str.substr(first_paren + 1, first_close - first_paren - 1);
  result.d1_sqr = std::stoul(d1_str);

  size_t second_start = first_close + 1;
  size_t second_paren = pair_str.find('(', second_start);
  if (second_paren == string::npos) {
    throw(std::invalid_argument(
        "Invalid NonInteractingPair format - missing second parentheses"));
  }

  string particle2_name =
      pair_str.substr(second_start, second_paren - second_start);

  size_t second_close = pair_str.find(')', second_paren);
  if (second_close == string::npos) {
    throw(std::invalid_argument("Invalid NonInteractingPair format - missing "
                                "second closing parentheses"));
  }

  string d2_str =
      pair_str.substr(second_paren + 1, second_close - second_paren - 1);
  result.d2_sqr = std::stoul(d2_str);

  // Find the decay channel that matches this particle pair
  bool found = false;
  for (uint dc = 0; dc < decay_channels.size(); ++dc) {
    const DecayChannelInfo& dci = decay_channels[dc];

    // Check if this decay channel matches the particle pair
    if ((dci.getParticle1Name() == particle1_name &&
         dci.getParticle2Name() == particle2_name) ||
        (dci.getParticle1Name() == particle2_name &&
         dci.getParticle2Name() == particle1_name)) {
      result.decay_channel_idx = dc;
      found = true;
      break;
    }
  }

  if (!found) {
    throw(std::invalid_argument("Particle pair '" + particle1_name + "," +
                                particle2_name +
                                "' not found in decay channels"));
  }

  return result;
}
