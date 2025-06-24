#include "chisq_detres.h"
#include "task_utils.h"
#include <filesystem>
#include <map>
#include <string>
using namespace std;

// *************************************************************

DeterminantResidualFit::DeterminantResidualFit(XMLHandler& xmlin,
                                               KBObsHandler* kboh,
                                               XMLHandler& xmlout,
                                               const string& outfile_stub)
    : ChiSquare() {
  KBOH = kboh;
  Kmat = 0;
  Kinv = 0;
  try {
    XMLHandler xmlf(xmlin, "DeterminantResidualFit");

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

    // speecify the folder name for the output
    // TODO: change outputs
    string outstub;
    xmlreadif(xmlin, "OutputStub", outstub, "DeterminantResidualFit");

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
    xmlreadif(xmlin, "OmegaMu", omega_mu, "DeterminantResidualFit");
    logger << "Omega mu = " << omega_mu << endl;

    //  get format for energies/mass (in terms of product with
    //  lattice time spacing or a ratio with reference mass)
    bool energy_ratios = true;
    {
      string reply;
      xmlreadif(xmlin, "DefaultEnergyFormat", reply, "DeterminantResidualFit");
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

    vector<KBObsInfo> labenergies;
    vector<MCEnsembleInfo> blockens;
    map<MCEnsembleInfo, uint> ensemble_idmap;

    list<XMLHandler> xmlkb = xmlf.find("KBBlock");
    uint blockcount = 0;
    uint numresiduals = 0;
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
        labenergies.push_back(KBObsInfo(mcens, rkey));
        nres++;
        numresiduals++;
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

    // get QuantizationCondition
    string qctype;
    xmlreadif(xmlf, "QuantizationCondition", qctype, "DeterminantResidualFit");
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

    //  connect files for input

    KBOH->connectSamplingFiles(sampfiles, needed_keys, verbose);
    logger << KBOH->getCurrentLog().str();
    KBOH->clearLog();

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

    if (numresiduals != labenergies.size())
      throw(std::runtime_error("Mismatch in residual sizes"));
    if (nsamplings == 0)
      throw(std::runtime_error(
          "Samplings numbers do not match for all ensembles"));

    Ecm_over_mref.resize(numresiduals);
    ensemble_id.resize(numresiduals);
    Bmat.resize(numresiduals);
    vector<const RVector*> labenergyvalues(numresiduals);
    vector<vector<RVector>> qcmsq_over_mrefsq(numresiduals);

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

      KBObsInfo lengthkey(mcens, MCObsInfo("LengthReference"));
      const RVector& mrefL = KBOH->getFullAndSamplingValues(lengthkey);
      KBObsInfo scalekey(mcens, MCObsInfo("KBScale"));
      const RVector& atrefmass = KBOH->getFullAndSamplingValues(scalekey);

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

      //  now evaluate Ecm/mref and the box matrix Bmat for each lab energy
      //  and each resampling  (form ratios if energy_ratio false)

      for (uint resind = indstart; resind < indstop; ++resind) {
        Ecm_over_mref[resind].resize(nsamp + 1);
        ensemble_id[resind] = mcensid;
        Bmat[resind].resize(nsamp + 1);
        qcmsq_over_mrefsq[resind].resize(nsamp + 1);
        const RVector* labenergy =
            &(KBOH->getFullAndSamplingValues(labenergies[resind]));
        if (!energy_ratios) {
          RVector labenergyratio(*labenergy);
          if (labenergyratio.size() != atrefmass.size())
            throw(std::runtime_error(
                "Size mismatch while forming energy ratios"));
          for (uint kk = 0; kk < labenergyratio.size(); ++kk)
            labenergyratio[kk] /= atrefmass[kk];
          labenergy = &(KBOH->putFullAndSamplings(labenergies[resind],
                                                  labenergyratio, true));
        }
        labenergyvalues[resind] = labenergy;
      }
      for (uint b = 0; b <= nsamp; ++b) {
        bqptr->setRefMassL(mrefL[b]);
        for (uint ci = 0; ci < bqptr->getNumberOfDecayChannels(); ++ci) {
          bqptr->setMassesOverRef(ci, (*(particlemass1[ci]))[b],
                                  (*(particlemass2[ci]))[b]);
        }
        for (uint resind = indstart; resind < indstop; ++resind) {
          const double& Elab = (*(labenergyvalues[resind]))[b];
          Ecm_over_mref[resind][b] = bqptr->getEcmOverMrefFromElab(Elab);
          bqptr->getQcmsqOverMrefsqFromElab(Elab, qcmsq_over_mrefsq[resind][b]);
          bqptr->getBoxMatrixFromEcm(Ecm_over_mref[resind][b], Bmat[resind][b]);
        }
      }
      for (uint resind = indstart; resind < indstop; ++resind) {
        logger << "   LabEnergyValue = " << (*(labenergyvalues[resind]))[0]
               << "  Ecm/mref = " << Ecm_over_mref[resind][0] << endl;
      }
      indstart += nres;
    }

    // output Ecm_over_mref, qcmsq_over_mref, Bmat samplings to file, if
    // requested
    if (!outfile_stub.empty()) {
      logger << "Outputting Ecm/mref, qcmsq/mrefsq, Box matrix elements to "
                "samplings files with stub "
             << outfile_stub << endl;
      std::vector<SamplingsPutHandler*> sampput(ensemble_idmap.size(), 0);
      bool overwrite = true;
      for (map<MCEnsembleInfo, uint>::const_iterator it =
               ensemble_idmap.begin();
           it != ensemble_idmap.end(); ++it) {
        string fname = outfile_stub;
        fname += ".hdf5[/ens" + make_string(it->second) + "]";
        SamplingsPutHandler* sp =
            new SamplingsPutHandler(KBOH->getBinsInfo(it->first),
                                    KBOH->getSamplingInfo(), fname, overwrite);
        sampput[it->second] = sp;
      }
      for (uint resind = 0; resind < numresiduals; ++resind) {
        uint ensid = ensemble_id[resind];
        const MCObsInfo& mkey = labenergies[resind].getMCObsInfo();
        // write Ecm/mref
        MCObsInfo Ecmkey(string("Ecm_mref[") + mkey.getObsName() + "]",
                         mkey.getObsIndex());
        sampput[ensid]->putData(Ecmkey, Ecm_over_mref[resind]);
        // write qcmsq/mrefsq for each channel
        uint nchan = qcmsq_over_mrefsq[resind][0].size();
        for (uint chan = 0; chan < nchan; ++chan) {
          MCObsInfo qcmsqkey(string("qcmsq_mrefsq[") + make_string(chan) +
                                 "][" + mkey.getObsName() + "]",
                             mkey.getObsIndex());
          uint nb = qcmsq_over_mrefsq[resind].size();
          RVector buff(nb);
          for (uint b = 0; b < nb; ++b)
            buff[b] = qcmsq_over_mrefsq[resind][b][chan];
          sampput[ensid]->putData(qcmsqkey, buff);
        }
      }
      // write out the Box matrix elements for each Elab
      uint indstart = 0;
      for (uint blocknum = 0; blocknum < nres_per_block.size(); ++blocknum) {
        uint nres = nres_per_block[blocknum];
        uint indstop = indstart + nres;
        BoxQuantization* bqptr = BQ[blocknum];
        for (uint row = 0; row < bqptr->getBasisSize(); ++row)
          for (uint col = 0; col < bqptr->getBasisSize(); ++col) {
            string belemname = bqptr->getKeyString(row, col);
            for (uint resind = indstart; resind < indstop; ++resind) {
              const MCObsInfo& ekey = labenergies[resind].getMCObsInfo();
              // TODO: revert quick/dirty fix for too long of names; just enumerate the energies
              MCObsInfo bkey(belemname + "[" + to_string(resind) + "]",
                             ekey.getObsIndex());
              uint nb = Bmat[resind].size();
              RVector buffre(nb), buffim(nb);
              for (uint b = 0; b < nb; ++b) {
                buffre[b] = Bmat[resind][b](row, col).real();
                buffim[b] = Bmat[resind][b](row, col).imag();
              }
              sampput[ensemble_id[resind]]->putData(bkey, buffre);
              bkey.setToImaginaryPart();
              sampput[ensemble_id[resind]]->putData(bkey, buffim);
            }
          }
        indstart += nres;
      }
      // close files
      for (uint k = 0; k < ensemble_idmap.size(); ++k) {
        delete sampput[k];
      }
    }

    KBOH->clearSamplings();
    initialize_base(numfitparams, numresiduals, nsamplings);
    xmlformat("DeterminantResidualConstruction", logger.str(), xmlout);
  } catch (const std::exception& xp) {
    string msg("Construction of DeterminantResidualFit failed: ");
    msg += xp.what();
    cout << msg << endl;
    clear();
    throw(std::runtime_error(msg));
  }
}

DeterminantResidualFit::~DeterminantResidualFit() { clear(); }

void DeterminantResidualFit::clear() {
  for (uint k = 0; k < BQ.size(); ++k)
    delete BQ[k];
  delete Kmat;
  delete Kinv;
  BQ.clear();
}

void DeterminantResidualFit::evalResidualsAndInvCovCholesky(
    const vector<double>& fitparams) {
  if (Kmat != 0)
    Kmat->setParameterValues(fitparams);
  else
    Kinv->setParameterValues(fitparams);
  RealSymmetricMatrix cov(nresiduals, 0.0);
  uint indstart = 0;
  bool bootstrap = KBOH->isBootstrapMode();
  vector<vector<double>> res(nresiduals,
                             vector<double>(Ecm_over_mref[0].size()));
  for (uint blocknum = 0; blocknum < nres_per_block.size(); ++blocknum) {
    uint nbres = nres_per_block[blocknum];
    BoxQuantization* bqptr = BQ[blocknum];
    uint nssize = Ecm_over_mref[indstart].size();
    for (uint k = 0; k < nbres; ++k) {
      uint kk = indstart + k;
      res[kk].resize(nssize);
      for (uint b = 0; b < nssize; ++b) {
        res[kk][b] = bqptr->getOmegaFromEcm(omega_mu, Ecm_over_mref[kk][b],
                                            Bmat[kk][b], qctype_enum).real();
      }
    }
    for (uint k = 0; k < nbres; ++k) {
      // cout << "Ecm_over_mref["<<k+indstart<<"] =
      // "<<Ecm_over_mref[k+indstart][resampling_index]<<endl;
      residuals[k + indstart] = res[k + indstart][resampling_index];
    }
    indstart += nbres;
  }

  // evaluate covariances; set covariances between different ensembles to zero
  for (uint k = 0; k < nresiduals; ++k)
    for (uint j = 0; j <= k; ++j) {
      if ((j == k) || (ensemble_id[k] == ensemble_id[j]))
        cov(k, j) = (bootstrap) ? KBOH->boot_covariance(res[k], res[j])
                                : KBOH->jack_covariance(res[k], res[j]);
    }

  CholeskyDecomposer CD;

  CD.getCholeskyOfInverse(cov, inv_cov_cholesky);
}

void DeterminantResidualFit::guessInitialFitParamValues(
    vector<double>& fitparams) const {
  if (Kmat != 0)
    fitparams = Kmat->getParameterValues();
  else
    fitparams = Kinv->getParameterValues();
}

void DeterminantResidualFit::getFitParamMCObsInfo(
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

void DeterminantResidualFit::do_output(XMLHandler& xmlout) const {
  xmlout.set_root("DeterminantResidualFit");
  XMLHandler xmlK;
  if (Kmat != 0)
    Kmat->output(xmlK);
  else
    Kinv->output(xmlK);
  xmlout.put_child(xmlK);
}

const std::vector<KFitParamInfo>&
DeterminantResidualFit::getFitParamInfos() const {
  if (Kmat != 0)
    return Kmat->getFitParameterInfos();
  else
    return Kinv->getFitParameterInfos();
}

// *********************************************************************
