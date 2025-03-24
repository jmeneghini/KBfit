#include "mcobs_info.h"
#include "xml_handler.h"
#include <ctime>
#include <map>

using namespace std;

void run2_an_obs(const MCObsInfo& mcobs) {
  cout << endl << endl << "  ***************************** " << endl << endl;
  cout << "Long output:" << endl;
  cout << mcobs.output(true) << endl;
  cout << "Short output:" << endl;
  string mcstr = mcobs.output(false);
  cout << mcstr << endl;

  cout << " isVacuum(): " << mcobs.isVacuum() << endl;
  cout << " isRealPart: " << mcobs.isRealPart() << endl;
  cout << " isImaginaryPart: " << mcobs.isImaginaryPart() << endl;
  cout << " isSimple(): " << mcobs.isSimple() << endl;
  cout << " isNonSimple(): " << mcobs.isNonSimple() << endl;
  cout << " isPrimary(): " << mcobs.isPrimary() << endl;
  cout << " isSecondary(): " << mcobs.isSecondary() << endl;

  try {
    string obsname(mcobs.getObsName());
    cout << "getObsName:" << obsname << endl;
  } catch (const std::exception& xx) {
    cout << "getObsName caught exception" << endl;
  }

  try {
    uint index = mcobs.getObsIndex();
    cout << "getObsIndex:" << index << endl;
  } catch (const std::exception& xx) {
    cout << "getObsIndex caught exception" << endl;
  }
}

void run_an_obs(const MCObsInfo& mcobs, const MCObsInfo& comparemcobs) {
  run2_an_obs(mcobs);
  cout << "mcobs1==mcobs2? :" << (mcobs == comparemcobs) << endl;
  cout << "mcobs1!=mcobs2? :" << (mcobs != comparemcobs) << endl;
  cout << "mcobs1<mcobs2? :" << (mcobs < comparemcobs) << endl;
  cout << endl << endl;
}

bool check_an_obs(XMLHandler& xmlobs) {
  cout << endl << endl << "  ***************************** " << endl << endl;
  cout << "Input XML: " << xmlobs.output() << endl << endl;
  try {
    MCObsInfo mcobs(xmlobs);
    cout << "Long output:" << endl;
    cout << mcobs.output(true) << endl;
    cout << "Short output:" << endl;
    string mcstr = mcobs.output(false);
    cout << mcstr << endl;

    bool vac = (xmlobs.count("Vacuum") == 1);
    bool secondary = (xmlobs.count("ObsName") == 1);
    bool primary = !secondary;
    bool simple = (vac) || (secondary && xmlobs.count("Simple") == 1);
    bool nonsimple = !simple;
    bool realpart;
    if (xmlobs.count("Arg") == 0)
      realpart = true;
    else {
      XMLHandler xmlarg(xmlobs, "Arg");
      string nv = xmlarg.get_text_content();
      realpart = (nv == "RealPart") || (nv == "Re");
    }
    bool imagpart = !realpart;

    bool result = true;
    bool checker;
    checker = mcobs.isVacuum();
    if (checker != vac) {
      cout << "vac mismatch" << endl;
      result = false;
    }
    checker = mcobs.isRealPart();
    if (checker != realpart) {
      cout << "realpart mismatch" << endl;
      result = false;
    }
    checker = mcobs.isImaginaryPart();
    if (checker != imagpart) {
      cout << "imagpart mismatch" << endl;
      result = false;
    }
    checker = mcobs.isSimple();
    if (checker != simple) {
      cout << "simple mismatch" << endl;
      result = false;
    }
    checker = mcobs.isNonSimple();
    if (checker != nonsimple) {
      cout << "nonsimple mismatch" << endl;
      result = false;
    }
    checker = mcobs.isPrimary();
    if (checker != primary) {
      cout << "primary mismatch" << endl;
      result = false;
    }
    checker = mcobs.isSecondary();
    if (checker != secondary) {
      cout << "secondary mismatch" << endl;
      result = false;
    }
    ComplexArg arg = (realpart) ? RealPart : ImaginaryPart;

    try {
      string obsname(mcobs.getObsName());
      if (!secondary) {
        cout << "problem getting obsname" << endl;
        result = false;
      }
    } catch (const std::exception& xx) {
      if (secondary) {
        cout << "problem getting obsname" << endl;
        result = false;
      }
    }

    try {
      uint index = mcobs.getObsIndex();
      if (!secondary) {
        cout << "problem getting obs index" << endl;
        result = false;
      }
      uint tcheck;
      xmlread(xmlobs, "Index", tcheck, "tester");
      if (index != tcheck) {
        cout << "problem getting obs index" << endl;
        result = false;
      }
    } catch (const std::exception& xx) {
      if (secondary) {
        cout << "problem getting obs index" << endl;
        result = false;
      }
    }

    try {
      string obsname(mcobs.getObsName());
      uint index = mcobs.getObsIndex();
      if (!secondary) {
        cout << "problem getting obsname" << endl;
        result = false;
      }
      MCObsInfo mcheck(obsname, index, simple, arg);
      if (mcheck != mcobs) {
        cout << "problem getting obsname" << endl;
        result = false;
      }
    } catch (const std::exception& xx) {
      if (secondary) {
        cout << "problem getting obsname" << endl;
        result = false;
      }
    }

    MCObsInfo mcheck2(mcobs);
    if (mcheck2 != mcobs) {
      cout << "constructor failed" << endl;
      result = false;
    }
    if (!vac) {
      if (realpart)
        mcheck2.setToImaginaryPart();
      else
        mcheck2.setToRealPart();
      if (mcheck2 == mcobs) {
        cout << "deliberate mismatch" << endl;
        result = false;
      }
    }
    mcheck2 = mcobs;
    if (mcheck2 != mcobs) {
      cout << "equality failed" << endl;
      result = false;
    }
    if (secondary) {
      mcheck2.resetObsIndex(mcobs.getObsIndex() + 8);
      if (mcheck2 == mcobs) {
        cout << "deliberate mismatch" << endl;
        result = false;
      }
    }
    mcheck2 = mcobs;
    if (mcheck2 != mcobs) {
      cout << "equality failed" << endl;
      result = false;
    }

    // cout << "numints = "<<mcobs.numints()<<endl;
    // cout << "number bytes = "<<mcobs.numbytes()<<endl;

    // unsigned int buf[24];
    // mcobs.copyTo(buf);

    // MCObsInfo mcobsBB(buf);
    // cout << "Check copyTo and create: "<<(mcobs==mcobsBB)<<endl;

    if (result)
      cout << "Test OK" << endl << endl;
    else
      cout << "Test FAILURE" << endl << endl;
    return result;
  } catch (const exception& xx) {
    cout << "Exception in creation: test FAILURE" << endl;
    cout << xx.what() << endl;
    return false;
  }
}

void testMCObsInfo(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMCObsInfo") == 0)
    return;

  cout << endl
       << endl
       << "***************************************************" << endl
       << endl;
  cout << "Testing MCObsInfo" << endl;

  MCObsInfo mc1;
  cout << mc1.output() << endl;

  list<XMLHandler> mcobsxml = xml_in.find("MCObservable");
  cout << "Found " << mcobsxml.size() << " MCObsInfo XML tags" << endl << endl;

  for (list<XMLHandler>::iterator it = mcobsxml.begin(); it != mcobsxml.end();
       it++) {
    try {
      cout << endl << endl << " ********************" << endl << endl;
      cout << "Input XML: " << it->output() << endl;
      MCObsInfo mc2(*it);
      cout << mc2.output(true) << endl;
      cout << mc2.output(false) << endl;
      cout << "Set to RealPart" << endl;
      mc2.setToRealPart();
      cout << mc2.output(false) << endl;
      cout << "Set to ImaginaryPart" << endl;
      mc2.setToImaginaryPart();
      cout << mc2.output(false) << endl;
    } catch (const std::invalid_argument& errmsg) {
      cout << "Whoops! Invalid XML: " << errmsg.what() << endl;
    }
  }

  mcobsxml = xml_in.find("MCObs");
  cout << "Found " << mcobsxml.size() << " MCObs XML tags" << endl << endl;

  for (list<XMLHandler>::iterator it = mcobsxml.begin(); it != mcobsxml.end();
       it++) {
    try {
      cout << endl << endl << " ********************" << endl << endl;
      cout << "Input XML: " << it->output() << endl;
      MCObsInfo mc2(*it);
      cout << mc2.output(true) << endl;
      cout << mc2.output(false) << endl;
      cout << "Set to RealPart" << endl;
      mc2.setToRealPart();
      cout << mc2.output(false) << endl;
      cout << "Set to ImaginaryPart" << endl;
      mc2.setToImaginaryPart();
      cout << mc2.output(false) << endl;
    } catch (const std::invalid_argument& errmsg) {
      cout << "Whoops! Invalid XML: " << errmsg.what() << endl;
    }
  }

  cout << endl
       << endl
       << " *************run an obs(mc1,mcs1)*****" << endl
       << endl;
  run_an_obs(mc1, mc1);

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

  uint ncheck = 0;
  uint nfail = 0;
  XMLHandler xmlcheck;
  for (uint simpcnt = 0; simpcnt < 2; simpcnt++)
    for (uint realcnt = 0; realcnt < 2; realcnt++)
      for (list<string>::iterator itn = obsnames.begin(); itn != obsnames.end();
           itn++) {
        for (list<uint>::iterator it = indices.begin(); it != indices.end();
             it++) {
          xmlcheck.set_root("MCObservable");
          xmlcheck.put_child("ObsName", *itn);
          xmlcheck.put_child("Index", make_string(*it));
          if (simpcnt > 0)
            xmlcheck.put_child("Simple");
          if (realcnt == 0)
            xmlcheck.put_child("Arg", "RealPart");
          else
            xmlcheck.put_child("Arg", "Im");
          if (!check_an_obs(xmlcheck))
            nfail++;
          ncheck++;
        }
      }

  cout << endl
       << endl
       << "Total number of checks = " << ncheck << endl
       << "Total number of check failures = " << nfail << endl
       << endl;
}
