#include "chisq_spectrum.h"
using namespace std;

SpectrumFit::SpectrumFit(XMLHandler& xmlin,
                          KBObsHandler* kboh,
                          XMLHandler& xmlout,
                          const string& outfile_stub) {
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
      // KBOH->setSamplingInfo(sampinfo);
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

    //  first set up the K or K^(-1) matrix // TODO: allow use of K or Kinv for cayley transfored QCs
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
    for (vector<DecayChannelInfo>::const_iterator dct = dcptr->begin();
         dct != dcptr->end(); dct++) {
      pnames.insert(dct->getParticle1Name());
      pnames.insert(dct->getParticle2Name());
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

    //  Loop over the KB quantization blocks to get the shift
    //  energies infos.

    vector<MCEnsembleInfo> blockens;
    map<MCEnsembleInfo, uint> ensemble_idmap;

    list<XMLHandler> xmlkb = xmlf.find("KBBlock");
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
      set<MCObsInfo>& kset = needed_keys[mcens];
      uint nres = 0;
      MCObsInfo rkey;
      list<XMLHandler> xmlee = it->find("EnergyShift");
      for (auto & eet : xmlee) {
        read_obs(eet, "EnergyShift", rkey, kset);
        KBObsInfo this_energy_key(mcens, rkey);
        energy_kobs_infos.emplace_back(this_energy_key);
        energy_samples_per_ensemble[ensemblecount].push_back(KBOH->getFullAndSamplingValues(this_energy_key));
        nres++;
      }
      if (nres == 0)
        throw(std::invalid_argument(
            "No energies available in at least one block"));
      n_energies_per_block.push_back(nres);
      blockcount++;
    }
    if (blockcount == 0) {
      throw(std::runtime_error("No data to analyze"));
    }

    // set the size of our vectors
    energy_samples_per_ensemble.resize(ensemblecount);
    mass_samples_per_ensemble.resize(ensemblecount);
    length_samples_per_ensemble.resize(ensemblecount);


    // get total number of energy residuals
    uint numenergies = energy_kobs_infos.size();
    //  connect files for input

    KBOH->connectSamplingFiles(sampfiles, needed_keys, verbose);
    logger << KBOH->getCurrentLog().str();
    KBOH->clearLog();

    // get QuantizationCondition
    string qctype;
    xmlreadif(xmlf, "QuantizationCondition", qctype, "DoPrint");
    if (qctype.empty()) {
      throw(std::invalid_argument("QuantizationCondition tag must be present"));
    }
    BoxQuantization* bqptr_dummy = BQ[0];
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

    // Bmat.resize(numresiduals);
    vector<const RVector*> labenergyvalues(numenergies);
    vector<vector<RVector>> qcmsq_over_mrefsq(numenergies);

    //  insert all of the fixed values into KBOH

    for (map<KBObsInfo, double>::const_iterator fx = fixed_values.begin();
         fx != fixed_values.end(); ++fx) {
      KBOH->putFixedValue(fx->first, fx->second, KBOH->getNumberOfResamplings(),
                          true);
    }

    //  read the reference mass time-spacing products for each ensemble into
    //  memory; evaluate reference lengths and particle masses for each ensemble
    uint current_ensemble_id = 0;
    for (map<MCEnsembleInfo, MCObsInfo>::iterator et = ref_at_mass.begin();
         et != ref_at_mass.end(); ++et) {
      const MCEnsembleInfo& mcens = et->first;
      uint nsamp = KBOH->getNumberOfResamplings();
      KBObsInfo atrefmasskey(mcens, et->second);
      const RVector& atrefmass0 = KBOH->getFullAndSamplingValues(atrefmasskey);
      if (atrefmass0.size() != (nsamp + 1))
        throw(std::runtime_error("Resampling size mismatch in KBfit"));
      KBObsInfo scalekey(mcens, MCObsInfo("KBScale"));
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
      MCObsInfo obskey("LengthReference", current_ensemble_id);
      KBObsInfo lengthkey(mcens, obskey);
      {
        RVector buff(nsamp + 1);
        if (anisotropy.find(mcens) == anisotropy.end()) { // isotropic case
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat);
        } else {
          KBObsInfo anisotropykey(mcens, anisotropy[mcens]);
          const RVector& anisotropy =
              KBOH->getFullAndSamplingValues(anisotropykey);
          if (anisotropy.size() != (nsamp + 1))
            throw(std::runtime_error("Resampling size mismatch in KBfit"));
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat) * anisotropy[k];
        }
        length_samples_per_ensemble.push_back(buff);
        KBOH->putFullAndSamplings(lengthkey, buff);
      }
      prior_obs_infos.push_back(obskey);
      prior_kobs_infos.push_back(lengthkey);

      //  get particle masses (form ratios if energy_ratio false)

      map<string, MCObsInfo>& pmap = particle_masses[mcens];
      for (set<string>::const_iterator pt = pnames.begin(); pt != pnames.end();
           ++pt) {
        MCObsInfo& current_particle_key = pmap[*pt];
        current_particle_key.resetObsIndex(current_ensemble_id);
        KBObsInfo masskey(mcens, current_particle_key);
        prior_obs_infos.push_back(current_particle_key);
        prior_kobs_infos.push_back(masskey);
        const RVector& mass = KBOH->getFullAndSamplingValues(
            masskey); // reads from file, gets into memory
        bool not_fixed = (fixed_values.find(masskey) == fixed_values.end());
        if ((!energy_ratios) && (not_fixed)) {
          RVector massratio(mass);
          if (massratio.size() != atrefmass.size())
            throw(
                std::runtime_error("Size mismatch while forming mass ratios"));
          for (uint kk = 0; kk < massratio.size(); ++kk)
            massratio[kk] /= atrefmass[kk];
          KBOH->putFullAndSamplings(masskey, massratio, true);
          mass_samples_per_ensemble[current_ensemble_id].push_back(massratio);
        } else {
          mass_samples_per_ensemble[current_ensemble_id].push_back(mass);
        }
      }
      current_ensemble_id++;
    }

    initialize_base(n_kmat_params, numenergies, nsamplings);
    
    // The SpectrumFit is designed to use a fixed covariance matrix calculated once.
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



SpectrumFit::~SpectrumFit() {
  clear();
}

void SpectrumFit::clear() {
  for (uint k = 0; k < BQ.size(); ++k)
    delete BQ[k];
  delete Kmat;
  delete Kinv;
  BQ.clear();
}

// Need to add non-Kmatrix fit parameters below here
// (L, mass, etc.), just use mean values.
void SpectrumFit::guessInitialFitParamValues(
    vector<double>& fitparams) const {
  // 1.  k-matrix parameters
  const std::vector<double>& k_fit_params =
      (Kmat ? Kmat->getParameterValues() : Kinv->getParameterValues());

  std::size_t offset = 0;                                  // where we're writing next
  std::copy(k_fit_params.begin(), k_fit_params.end(),      // copy K-matrix params
            fitparams.begin() + offset);
  offset += k_fit_params.size();                           // advance offset

  // 2.  Ensemble-specific parameters:  L, m01 … m0N,  L, m11 … etc.
  // We'll use their full values to start with.
  for (std::size_t e = 0; e < length_samples_per_ensemble.size(); ++e) {

    // L for ensemble e
    fitparams[offset++] = length_samples_per_ensemble[e][0];

    // masses for ensemble e
    const auto& masses = mass_samples_per_ensemble[e];
    for (const auto& mass_samples : masses) {
      fitparams[offset++] = mass_samples[0];
    }
  }
}

void SpectrumFit::getFitParamMCObsInfo(
    vector<MCObsInfo>& fitinfos) const {
  const vector<KFitParamInfo>& k_infos =
    (Kmat ? Kmat->getFitParameterInfos() : Kinv->getFitParameterInfos());
  uint offset = 0;
  for (const auto& info : k_infos)
    fitinfos[offset++] = MCObsInfo(info.getMCObsName());
  for (const auto& obs : prior_obs_infos)
    fitinfos[offset++] = obs;
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

// This method calculates residuals. Covariance is fixed after initializeInvCovCholesky.
void SpectrumFit::evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) {
  // Think im forced to copy for the kmat parameters, unless we're > c++ 20 (std::span)
  std::vector<double> k_params(fitparams.begin(), fitparams.begin() + n_kmat_params);
  (Kmat ? Kmat->setParameterValues(k_params) : Kinv->setParameterValues(k_params));
  std::vector<double> prior_params(fitparams.begin() + n_kmat_params, fitparams.end());

  uint offset = 0;
  uint block_offset = 0;
  for (uint ensemble = 0; ensemble < length_samples_per_ensemble.size(); ++ensemble) {
    // setup the BQ for the ensemble/block
    // first eval the length residual
    double length_param = prior_params[offset];
    residuals[offset]
        = length_samples_per_ensemble[ensemble][resampling_index] - length_param;
    offset++;
    uint mass_offset = 0;
    for (uint decay_channel = 0; decay_channel < n_decay_channels; ++decay_channel) {
      double mass1 = prior_params[offset];
      residuals[offset]
          = mass_samples_per_ensemble[ensemble][mass_offset++][resampling_index] - mass1;
      offset++;
      double mass2;
      if (are_decay_channels_identical[decay_channel]) {
        mass2 = mass1;
      }
      else {
        mass2 = prior_params[offset];
        residuals[offset]
            = mass_samples_per_ensemble[ensemble][mass_offset++][resampling_index] - mass2;
        offset++;
      }
      for (uint block = 0; block < n_blocks_per_ensemble[ensemble ]; ++block) {
        // set all the masses and lengths for all BQ blocks in the ensemble
        BQ[block + block_offset]->setRefMassL(length_param);
        BQ[block + block_offset]->setMassesOverRef(decay_channel, mass1, mass2);
      }
    }
    // now BQs are set up for the ensemble, so we can get the energy predictions
    uint energy_offset = 0;
    for (uint block = 0; block < n_blocks_per_ensemble[ensemble]; ++block) {
      BoxQuantization* this_block_bq = BQ[block + block_offset];
      std::vector<double> energy_shift_predictions;
      energy_shift_predictions.reserve(n_energies_per_block[block + block_offset]); // getRoots uses this.
      std::vector<uint> fn_calls;
      this_block_bq->getDeltaERootsInElabInterval(omega_mu, Elab_min, Elab_max,
                                      qctype_enum, root_finder_config,
                                      energy_shift_predictions, fn_calls);
      // both the energy obs and predictions are sorted by increasing Ecm
      for (uint energy_index = 0; energy_index < n_energies_per_block[block + block_offset]; ++energy_index) {
        residuals[offset++]
          = energy_samples_per_ensemble[ensemble][energy_offset++][resampling_index]
            - energy_shift_predictions[energy_index];
      }
    }
    block_offset += n_blocks_per_ensemble[ensemble];
  }
  // InvCovCholesky is already initialized and hasn't changed,
  // so we don't need to recompute it here.
}



// By introducing the priors for our MC observables,
// cov(r_i, r_j) simplifies to cov(R_i, R_j), where
// R is the observable
void SpectrumFit::initializeInvCovCholesky() {
  std::vector<const RVector*> obs_samples;
  std::vector<MCEnsembleInfo> obs_ensemble_infos;
  obs_samples.reserve(nresiduals);
  obs_ensemble_infos.reserve(nresiduals);

  uint block_offset = 0;

  for (uint ens = 0; ens < length_samples_per_ensemble.size(); ++ens) {
    // length
    obs_samples.push_back( &length_samples_per_ensemble[ens]            );
    obs_ensemble_infos.push_back( MCEnsembleInfo(ens) );  // or whatever ctor

    // masses
    uint mass_offset = 0;
    for (uint dc = 0; dc < n_decay_channels; ++dc)
    {
      obs_samples.push_back(
          &mass_samples_per_ensemble[ens][mass_offset++] );
      obs_ensemble_infos.push_back( MCEnsembleInfo(ens) );

      if (!are_decay_channels_identical[dc])
      {
        obs_samples.push_back(
            &mass_samples_per_ensemble[ens][mass_offset++] );
        obs_ensemble_infos.push_back( MCEnsembleInfo(ens) );
      }
    }

    // -- energy Delta E’s --------------------------------------------------------
    uint energy_offset = 0;
    for (uint block = 0; block < n_blocks_per_ensemble[ens]; ++block)
    {
      uint nE = n_energies_per_block[block + block_offset];
      for (uint e = 0; e < nE; ++e)
      {
        obs_samples.push_back(
            &energy_samples_per_ensemble[ens][energy_offset++] );
        obs_ensemble_infos.push_back( MCEnsembleInfo(ens) );
      }
    }
    block_offset += n_blocks_per_ensemble[ens];
  }

  RealSymmetricMatrix cov(nresiduals, 0.0);
  bool bootstrap_mode = KBOH->isBootstrapMode();

  for (uint k = 0; k < nresiduals; ++k) {
    for (uint j = 0; j <= k; ++j) {
      // Set covariance to zero if observables are from different ensembles
      if (obs_ensemble_infos[k] == obs_ensemble_infos[j]) {
        cov(k, j) = bootstrap_mode ? KBOH->boot_covariance(*(obs_samples[k]), *(obs_samples[j]))
                                   : KBOH->jack_covariance(*(obs_samples[k]), *(obs_samples[j]));
      }
      else {
        cov(k, j) = 0.0; // different ensembles, 0 covariance
      }
    }
  }

  CholeskyDecomposer CD;
  CD.getCholeskyOfInverse(cov, inv_cov_cholesky);
}

