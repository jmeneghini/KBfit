#include "kbobs_handler.h"
#include "xml_handler.h"
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <random>
using namespace std;

#ifdef SINGLEPRECISION
const double eps = 1e-5;
const double ceps = 5e-5;
#else
const double eps = 1e-11;
const double ceps = 5e-11;
#endif
const double deps = 1e-11;

double get_random_double() {
  return double(rand() % 1048576) / 1048576.0 - 0.5;
}

void testMCObsGetHandler(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMCObsGetHandler") == 0)
    return;

  cout << endl << "Starting TestMCObsGetHandler" << endl;

  try {
    XMLHandler xmlr(xml_in, "TestMCObsGetHandler");
    MCSamplingInfo sampinfo(xmlr);
    MCObsGetHandler MCOGH(xmlr);
    cout << "done MCOGH set up" << endl;

    cout << "is Jackknife mode? " << MCOGH.isJackknifeMode() << endl;
    cout << "is Bootstrap mode? " << MCOGH.isBootstrapMode() << endl;

    set<MCEnsembleInfo> ensset = MCOGH.getEnsembleInfos();
    for (set<MCEnsembleInfo>::const_iterator st = ensset.begin();
         st != ensset.end(); ++st) {
      cout << st->output() << endl;
      const MCBinsInfo& binfo = MCOGH.getBinsInfo(*st);
      const MCSamplingInfo& sampinfo = MCOGH.getSamplingInfo();
      cout << "BinInfo: " << binfo.output() << endl;
      cout << "SamplingInfo: " << sampinfo.output() << endl;
    }

    RVector samp;
    XMLHandler xmlg(xmlr, "GetTests");
    list<XMLHandler> kbxml = xmlg.find("KBObservable");
    for (list<XMLHandler>::iterator it = kbxml.begin(); it != kbxml.end();
         ++it) {
      KBObsInfo kbkey(*it);
      cout << "KBkey:" << kbkey.output();
      cout << "query = " << MCOGH.querySamplings(kbkey) << endl;
      try {
        MCOGH.getSamplings(kbkey, samp);
        cout << "samp size = " << samp.size() << "  samp[0] = " << samp[0]
             << endl;
      } catch (...) {
        cout << "read could not happen" << endl;
      }
      samp.clear();
      bool flag = MCOGH.getSamplingsMaybe(kbkey, samp);
      if (flag)
        cout << "samp size = " << samp.size() << "  samp[0] = " << samp[0]
             << endl;
      else
        cout << "could not read" << endl;
    }

    cout << endl
         << "************************************************" << endl
         << endl;
    cout << "Now testing KCObsHandler" << endl;
    KBObsHandler KBOH(sampinfo);
    cout << "get sampling mode = " << KBOH.getSamplingMode() << endl;
    cout << "sampling info = " << KBOH.getSamplingInfo().output() << endl;
    cout << "is jackknife mode? " << KBOH.isJackknifeMode() << endl;
    cout << "is bootstrap mode? " << KBOH.isBootstrapMode() << endl;

    for (list<XMLHandler>::iterator it = kbxml.begin(); it != kbxml.end();
         ++it) {
      KBObsInfo kbkey(*it);
      cout << "KBkey:" << kbkey.output();
      try {
        cout << "MCBinsInfo: " << KBOH.getBinsInfo(kbkey).output() << endl;
        cout << "number of resamplings = " << KBOH.getNumberOfResamplings()
             << endl;
      } catch (...) {
        cout << "no getBinsInfo" << endl;
      }
      cout << "query = " << KBOH.queryFullAndSamplings(kbkey) << endl;
      try {
        const RVector& res = KBOH.getFullAndSamplingValues(kbkey);
        cout << "res size = " << res.size() << "  res[0] = " << res[0] << endl;
      } catch (...) {
        cout << "read could not happen" << endl;
      }
      // KBOH.clearSamplings();
      const RVector* result = KBOH.getFullAndSamplingValuesMaybe(kbkey);
      if (result == 0)
        cout << "Maybe read do not happen" << endl;
      else
        cout << "result[4]=" << (*result)[4] << endl;
      try {
        cout << "full value = " << KBOH.getFullSampleValue(kbkey) << endl;
        cout << "index 4 value = " << KBOH.getSamplingValue(kbkey, 4) << endl;
      } catch (...) {
        cout << "get values failed" << endl;
      }
      double rr;
      if (KBOH.getSamplingValueMaybe(kbkey, 4, rr))
        cout << "rr = " << rr << endl;
      else
        cout << "rr could not be read" << endl;
    }

    cout << endl << "MCEnsembleInfos involved:" << endl;
    set<MCEnsembleInfo> kset(KBOH.getEnsembleInfos());
    for (set<MCEnsembleInfo>::const_iterator kt = kset.begin();
         kt != kset.end(); ++kt)
      cout << kt->output() << endl;

    KBObsInfo kbkey(MCEnsembleInfo("indep"), MCObsInfo("Kangaroo", 8));
    KBOH.putSamplingValue(kbkey, 0, 6.000, true);

  } catch (const std::exception& err) {
    cerr << "  Error: " << err.what() << endl;
    cerr << "  ... exiting..." << endl;
    exit(1);
  }
}

// ***********************************************************************

/*
bool compare_vectors(const RVector& v1, const RVector& v2, double factor)
{
 if (v1.size()!=v2.size()) return false;
 for (uint k=0;k<v1.size();k++)
    if (std::abs(v1[k]-v2[k])>factor*std::max(std::abs(0.5*(v1[k]+v2[k])),1.0))
       return false;
 return true;
}


bool compare_floats(const double& v1, const double& v2, double factor)
{
 if (std::abs(v1-v2)>factor*std::max(std::abs(0.5*(v1+v2)),1.0)){ cout
<<"MISMATCH: "<< v1<<" "<<v2<<endl; return false;} return true;
}



void samplings_data_assign(double& data, int serind, double Aval, double Cval)
{
 double g=Aval*(1.0+Cval*serind)*get_random_double();
 data=0.4172*g;
}


void make_fake_samplings_file(const MCObsInfo& obskey, const std::string&
filename, const MCBinsInfo& bininfo, const MCSamplingInfo& sampinfo, double
Aval, double Cval, SamplingsCorrect& SC)
{
 SamplingsPutHandler SP(bininfo,sampinfo,filename,true,false);
 double data;
 uint nsamp=sampinfo.getNumberOfReSamplings();
 vector<double> buffer;
   // first is full sampling, then do the samplings
 for (int serind=0;serind<=int(nsamp);++serind){
    samplings_data_assign(data,serind,Aval,Cval);
    buffer.push_back(data);}
 SP.putData(obskey,RVector(buffer));
 SC.addValue(obskey,RVector(buffer));
}


void SamplingsCorrect::addValue(const MCObsInfo& mcobs, const RVector& insamp)
{
 map<MCObsInfo,RVector >::iterator it=values.find(mcobs);
 if (it!=values.end())
    it->second=insamp;
 else{
    values.insert(make_pair(mcobs,insamp));}
}

bool SamplingsCorrect::getCorrect(const MCObsInfo& mcobs, RVector& result) const
{
 map<MCObsInfo,RVector >::const_iterator it=values.find(mcobs);
 if (it==values.end()){ result.clear(); return false;}
 result=it->second;
 return true;
}

*/

// *******************************************************************

void testMCObsGetHandlerFake(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMCObsGetHandlerFake") == 0)
    return;

  cout << endl << "Starting TestMCObsGetHandlerFake" << endl;
  srand(time(NULL));

  std::default_random_engine generator;

  try {
    XMLHandler xmlr(xml_in, "TestMCObsGetHandlerFake");
    list<XMLHandler> getfilexml = xmlr.find("MakeAFile");
    for (list<XMLHandler>::iterator it = getfilexml.begin();
         it != getfilexml.end(); ++it) {
      MCSamplingInfo sampinfo(*it);
      MCBinsInfo binfo(*it);
      string filename;
      xmlread(*it, "FileName", filename, "TestMCObsGetHandlerFake");
      SamplingsPutHandler SPH(binfo, sampinfo, filename, true, false);
      list<XMLHandler> obsxml = it->find("DataItem");
      for (list<XMLHandler>::iterator ot = obsxml.begin(); ot != obsxml.end();
           ++ot) {
        MCObsInfo obskey(*ot);
        unsigned int nsamp = sampinfo.getNumberOfReSamplings();
        double mean, stddev;
        xmlread(*ot, "Mean", mean, "TestFake");
        xmlread(*ot, "StdDev", stddev, "TestFake");
        Vector<double> data(nsamp + 1);
        std::normal_distribution<double> distribution(mean, stddev);
        for (uint k = 0; k <= nsamp; ++k)
          data[k] = distribution(generator);
        SPH.putData(obskey, data);
      }
    }
  } catch (const std::exception& err) {
    cerr << "  Error: " << err.what() << endl;
    cerr << "  ... exiting..." << endl;
    exit(1);
  }
}
