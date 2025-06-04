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

    //  first set up the K or K^(-1) matrix
    uint numfitparams;
    int k1 = xmlf.count_among_children("KtildeMatrix");
    int k2 = xmlf.count_among_children("KtildeMatrixInverse");
    if ((k1 + k2) != 1)
      throw(std::invalid_argument(
          "A single KtildeMatrix or KtildeMatrixInverse tag must be present"));
    const vector<DecayChannelInfo>* dcptr;
    if (k1 == 1) {
      Kmat = new KtildeMatrixCalculator(xmlf, true);
      numfitparams = Kmat->getNumberOfParameters();
      dcptr = &(Kmat->getDecayChannelInfos());
    } else {
      Kinv = new KtildeInverseCalculator(xmlf, true);
      numfitparams = Kinv->getNumberOfParameters();
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

    //  Loop over the KB quantization blocks to get the lab-frame
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
      list<XMLHandler> xmlee = it->find("LabFrameEnergy");
      for (list<XMLHandler>::iterator eet = xmlee.begin(); eet != xmlee.end();
           ++eet) {
        read_obs(*eet, "LabFrameEnergy", rkey, kset);
        energy_obs_infos.push_back(KBObsInfo(mcens, rkey));
        nres++;
      }
      if (nres == 0)
        throw(std::invalid_argument(
            "No energies available in at least one block"));
      nres_per_block.push_back(nres);
      blockcount++;
    }
    if (blockcount == 0) {
      throw(std::runtime_error("No data to analyze"));
    }

    // get total number of residuals
    uint numenergies = energy_obs_infos.size();
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
              KBOH->getFullAndSamplingValues(anisotropykey);
          if (anisotropy.size() != (nsamp + 1))
            throw(std::runtime_error("Resampling size mismatch in KBfit"));
          for (uint k = 0; k <= nsamp; ++k)
            buff[k] = atrefmass[k] * double(Llat) * anisotropy[k];
        }
        KBOH->putFullAndSamplings(lengthkey, buff);
      }

      //  get particle masses (form ratios if energy_ratio false)

      map<string, MCObsInfo>& pmap = particle_masses[mcens];
      for (set<string>::const_iterator pt = pnames.begin(); pt != pnames.end();
           ++pt) {
        KBObsInfo masskey(mcens, pmap[*pt]);
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
        }
      }
    }

    //  now evaluate Ecm and the box matrices
    //  for each block and each resampling

    uint indstart = 0;
    for (uint blocknum = 0; blocknum < nres_per_block.size(); ++blocknum) {
      uint nres = nres_per_block[blocknum];
      if (verbose) {
        logger << "Setting up box matrix and Ecm for block number " << blocknum
               << endl;
        logger << "Number of energy levels = " << nres << endl;
      }
      uint indstop = indstart + nres;
      BoxQuantization* bqptr = BQ[blocknum];
      logger << "MomRay " << bqptr->getMomRay()
             << "   P^2 = " << bqptr->getTotalMomentumIntegerSquared()
             << " Box Irrep " << bqptr->getLittleGroupBoxIrrep() << endl;
      const MCEnsembleInfo& mcens = blockens[blocknum];
      uint mcensid = ensemble_idmap[mcens];
      uint nsamp = KBOH->getNumberOfResamplings();
      map<string, MCObsInfo>& pmap = particle_masses[mcens];

      // put all priors into vector
      KBObsInfo lengthkey(mcens, MCObsInfo("LengthReference"));
      prior_obs_infos.push_back(lengthkey);


      uint nchan = bqptr->getNumberOfDecayChannels();
      vector<const RVector*> particlemass1(nchan), particlemass2(nchan);
      for (uint ci = 0; ci < nchan; ++ci) {
        const DecayChannelInfo& chan = bqptr->getDecayChannelInfo(ci);
        const string& pname1 = chan.getParticle1Name();
        const string& pname2 = chan.getParticle2Name();
        KBObsInfo mass1key(mcens, pmap[pname1]);
        KBObsInfo mass2key(mcens, pmap[pname2]);
        particlemass1[ci] = &(KBOH->getFullAndSamplingValues(mass1key));
        particlemass2[ci] = &(KBOH->getFullAndSamplingValues(mass2key));
      }
    }

    KBOH->clearSamplings();
    initialize_base(numfitparams, numenergies, nsamplings);
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
  if (Kmat != 0)
    fitparams = Kmat->getParameterValues();
  else
    fitparams = Kinv->getParameterValues();
}

void SpectrumFit::getFitParamMCObsInfo(
    vector<MCObsInfo>& fitinfos) const {
  const vector<KFitParamInfo>* fpptr = 0;
  if (Kmat != 0)
    fpptr = &(Kmat->getFitParameterInfos());
  else
    fpptr = &(Kinv->getFitParameterInfos());
  fitinfos.resize(fpptr->size());
  for (uint k = 0; k < fitinfos.size(); ++k) {
    fitinfos[k] = MCObsInfo((*fpptr)[k].getMCObsName());
  }
}

// might want to modify this at some point once I
// better understand its purpose
void SpectrumFit::do_output(XMLHandler& xmlout) const {
  xmlout.set_root("DeterminantResidualFit");
  XMLHandler xmlK;
  if (Kmat != 0)
    Kmat->output(xmlK);
  else
    Kinv->output(xmlK);
  xmlout.put_child(xmlK);
}

void SpectrumFit::evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams) {

}


// By introducing the priors for our MC observables,
// cov(r_i, r_j) simplifies to cov(R_i, R_j), where
// R is the observable
void SpectrumFit::initializeInvCov() {
  bool bootstrap = KBOH->isBootstrapMode();

  // get params
  vector<MCObsInfo> param_infos;
  getFitParamMCObsInfo(param_infos);
  vector<vector<double>> obs_samples(param_infos.size());

  RealSymmetricMatrix cov(nresiduals, 0.0);
  for (uint k = 0; k < nresiduals; ++k)
    for (uint j = 0; j <= k; ++j) {
      if ((j == k) || (Ecm_ensemble_id[k] == Ecm_ensemble_id[j]))
        cov(k, j) = (bootstrap) ? KBOH->boot_covariance(res[k], res[j])
                                : KBOH->jack_covariance(res[k], res[j]);
    }
}

void SpectrumFit::read_obs(XMLHandler& xmlin, const string& tag,
                                      bool get_name, MCObsInfo& obskey,
                                      set<MCObsInfo>& kset, string& name,
                                      const MCEnsembleInfo& mcens,
                                      map<KBObsInfo, double>& fixed_values) {
  try {
    XMLHandler xmlt(xmlin, tag);
    name.clear();
    if (get_name) {
      xmlread(xmlt, "Name", name, tag);
    }
    uint mcount = xmlt.count("MCObs") + xmlt.count("MCObservable");
    uint fcount = xmlt.count("FixedValue");
    if ((mcount + fcount) != 1)
      throw(std::invalid_argument("No MCObs/MCObservable or FixedValue"));
    if (mcount == 1) {
      obskey = MCObsInfo(xmlt);
      if ((obskey.isImaginaryPart()) || (obskey.isSimple()))
        throw(std::invalid_argument("MCObsInfo must be nonsimple and real"));
      if (obskey == MCObsInfo("KBScale"))
        throw(std::invalid_argument(
            "KBScale is reserved and cannot be an input MCObsInfo"));
      if (kset.find(obskey) != kset.end())
        throw(std::invalid_argument("duplicate MCObsInfo"));
      kset.insert(obskey);
    } else {
      double fixedvalue;
      xmlread(xmlt, "FixedValue", fixedvalue, tag);
      string kbname(tag);
      if (!name.empty())
        kbname += "_" + name;
      obskey = MCObsInfo(kbname);
      KBObsInfo kbkey(mcens, obskey);
      fixed_values.insert(make_pair(kbkey, fixedvalue));
    }
  } catch (const std::exception& xp) {
    string msg = string("For tag ") + tag;
    throw(std::invalid_argument(msg + string(": ") + xp.what()));
  }
}

//  This routine does the following:
//    -- searches only within an XML tag with name specified in "tag"
//    -- an <MCObs> or <MCObservable> must be encountered, then "obskey"
//        is assigned this tag; if "obskey" is already in "kset", an
//        exception is thrown, but if not, then "obskey" is inserted
//        into "kset".
//    -- "obskey" must be nonsimple and real.

void SpectrumFit::read_obs(XMLHandler& xmlin, const string& tag,
                                      MCObsInfo& obskey, set<MCObsInfo>& kset) {
  try {
    XMLHandler xmlt(xmlin, tag);
    uint mcount = xmlt.count("MCObs") + xmlt.count("MCObservable");
    uint fcount = xmlt.count("FixedValue");
    if ((mcount != 1) || (fcount != 0))
      throw(std::invalid_argument(
          "No MCObs/MCObservable or disallowed FixedValue"));
    obskey = MCObsInfo(xmlt);
    if ((obskey.isImaginaryPart()) || (obskey.isSimple()))
      throw(std::invalid_argument("MCObsInfo must be nonsimple and real"));
    if (obskey == MCObsInfo("KBScale"))
      throw(std::invalid_argument(
          "KBScale is reserved and cannot be an input MCObsInfo"));
    if (kset.find(obskey) != kset.end())
      throw(std::invalid_argument("duplicate MCObsInfo"));
    kset.insert(obskey);
  } catch (const std::exception& xp) {
    string msg = string("For tag ") + tag;
    throw(std::invalid_argument(msg + string(": ") + xp.what()));
  }
}

