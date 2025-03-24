#include "kbobs_info.h"
#include "xml_handler.h"
#include <ctime>
#include <map>

using namespace std;

void run2_an_obs(const KBObsInfo& kbobs) {
  cout << endl << endl << "  ***************************** " << endl << endl;
  cout << "Long output:" << endl;
  cout << kbobs.output(true) << endl;
  cout << "Short output:" << endl;
  string kbstr = kbobs.output(false);
  cout << kbstr << endl;
}

void run_an_obs(const KBObsInfo& kbobs, const KBObsInfo& comparekbobs) {
  run2_an_obs(kbobs);
  cout << "kbobs1==kbobs2? :" << (kbobs == comparekbobs) << endl;
  cout << "kbobs1!=kbobs2? :" << (kbobs != comparekbobs) << endl;
  cout << "kbobs1<kbobs2? :" << (kbobs < comparekbobs) << endl;
  cout << endl << endl;
}

void testKBObsInfo(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestKBObsInfo") == 0)
    return;

  cout << endl
       << endl
       << "***************************************************" << endl
       << endl;
  cout << "Testing KBObsInfo" << endl;

  list<XMLHandler> kbobsxml = xml_in.find("KBObservable");
  cout << "Found " << kbobsxml.size() << " KBObservable XML tags" << endl
       << endl;

  MCEnsembleInfo ens0("clover_s16_t128_ud840_s743");

  for (list<XMLHandler>::iterator it = kbobsxml.begin(); it != kbobsxml.end();
       it++) {
    try {
      cout << endl << endl << " ********************" << endl << endl;
      cout << "Input XML: " << it->output() << endl;
      KBObsInfo kb2(*it);
      cout << kb2.output(true) << endl;
      cout << kb2.output(false) << endl;
      cout << "Set to RealPart" << endl;
      kb2.setToRealPart();
      cout << kb2.output(false) << endl;
      cout << "Set to ImaginaryPart" << endl;
      kb2.setToImaginaryPart();
      cout << kb2.output(false) << endl;
      cout << "Reset index to 763" << endl;
      kb2.resetObsIndex(763);
      cout << kb2.output(false) << endl;
      cout << "Reset ensemble to clover_s16_t128_ud840_s743" << endl;
      kb2.resetMCEnsembleInfo(ens0);
      cout << kb2.output(false) << endl;
    } catch (const std::invalid_argument& errmsg) {
      cout << "Whoops! Invalid XML: " << errmsg.what() << endl;
    }
  }
  /*
   cout << endl<<endl<<" *************run an obs(mc1,mcs1)*****"<<endl<<endl;
   run_an_obs(mc1,mc1);

   list<string> obsnames;
   obsnames.push_back("Kincardine");
   obsnames.push_back("Goderich");
   obsnames.push_back("Toronto834%!#");
   obsnames.push_back("NiagaraFalls");
   list<uint> indices;
   indices.push_back(4);
   indices.push_back(22);
   indices.push_back(11);
   indices.push_back(39);
   indices.push_back(2);
   indices.push_back(8);
   indices.push_back(1);

   uint ncheck=0;
   uint nfail=0;
   XMLHandler xmlcheck;
   for (uint simpcnt=0;simpcnt<2;simpcnt++)
   for (uint realcnt=0;realcnt<2;realcnt++)
   for (list<string>::iterator itn=obsnames.begin();itn!=obsnames.end();itn++){
   for (list<uint>::iterator it=indices.begin();it!=indices.end();it++){
      xmlcheck.set_root("MCObservable");
      xmlcheck.put_child("ObsName",*itn);
      xmlcheck.put_child("Index",make_string(*it));
      if (simpcnt>0) xmlcheck.put_child("Simple");
      if (realcnt==0) xmlcheck.put_child("Arg","RealPart");
      else xmlcheck.put_child("Arg","Im");
      if (!check_an_obs(xmlcheck)) nfail++;
      ncheck++;}}


   cout << endl<<endl<<"Total number of checks = "<<ncheck<<endl
        <<"Total number of check failures = "<<nfail<<endl<<endl;

  */
}
