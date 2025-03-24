#include "K_matrix_info.h"
#include <fstream>
using namespace std;

bool invalidJLS(uint Jtimestwo, uint L, uint Stimestwo) {
  return (((Jtimestwo % 2) != (Stimestwo % 2)) ||
          ((2 * L) < ((Jtimestwo >= Stimestwo) ? (Jtimestwo - Stimestwo)
                                               : (Stimestwo - Jtimestwo))) ||
          ((2 * L) > (Jtimestwo + Stimestwo)));
}

bool invalidJLS2(uint Jtimestwo, uint L, uint Stimestwo) {
  uint L2 = 2 * L;
  uint J2min = (L2 >= Stimestwo) ? (L2 - Stimestwo) : (Stimestwo - L2);
  uint J2max = L2 + Stimestwo;
  for (uint Jg = J2min; Jg <= J2max; Jg += 2)
    if (Jtimestwo == Jg)
      return false;
  return true;
}

void do_Kfitparamtests(XMLHandler& xmlin) {

  for (uint L = 0; L < 16; L++)
    for (uint S2 = 0; S2 < 16; S2++)
      for (uint a = 0; a < 16; a++) {
        KIndex K1(L, S2, a);
        cout << K1.str() << endl;
        uint LL = K1.getL();
        uint SS = K1.getStimestwo();
        uint aa = K1.getChannelIndex();
        if ((LL != L) || (S2 != SS) || (a != aa))
          cout << "ERROR" << endl;
      }

  cout << endl << endl << "Testing invalidJLS" << endl;
  for (uint J2 = 0; J2 < 9; J2++)
    for (uint L = 0; L < 16; L++)
      for (uint S2 = 0; S2 < 8; S2++) {
        bool chk = invalidJLS(J2, L, S2);
        bool chk2 = invalidJLS2(J2, L, S2);
        if (!chk)
          cout << "JLS " << double(J2) / 2.0 << " " << L << " "
               << double(S2) / 2.0 << endl;
        if (chk != chk2)
          cout << "ERROR in invalidJLS" << endl;
      }
  cout << "Test done" << endl;

  cout << "Testing KElementInfo" << endl << endl;

  uint imax = 7;
  uint count1 = 0, count2 = 0;
  for (uint J2 = 0; J2 < imax; J2++)
    for (uint Lp = 0; Lp < imax; Lp++)
      for (uint S2p = 0; S2p < imax; S2p++)
        for (uint ap = 0; ap < 4; ap++)
          for (uint L = 0; L < imax; L++)
            for (uint S2 = 0; S2 < imax; S2++)
              for (uint a = 0; a < 4; a++) {
                try {
                  KElementInfo K(J2, Lp, S2p, ap, L, S2, a);
                  KElementInfo K2(J2, L, S2, a, Lp, S2p, ap);
                  if (K != K2)
                    cout << "ERROR" << endl;
                  if (K.getJtimestwo() != J2)
                    cout << "J2 error" << endl;
                  KIndex row(K.getRow());
                  KIndex inrow(Lp, S2p, ap);
                  KIndex col = K.getColumn();
                  KIndex incol(L, S2, a);
                  if ((row == inrow) && (col == incol)) {
                    if ((K.getRowL() != Lp) || (K.getRowStimestwo() != S2p) ||
                        (K.getRowChannelIndex() != ap) ||
                        (K.getColumnL() != L) ||
                        (K.getColumnStimestwo() != S2) ||
                        (K.getColumnChannelIndex() != a))
                      cout << "ERROR in indices: " << endl;
                    if ((row.getL() != Lp) || (row.getStimestwo() != S2p) ||
                        (row.getChannelIndex() != ap))
                      cout << "ERROR in row indices: " << endl;
                    if ((col.getL() != L) || (col.getStimestwo() != S2) ||
                        (col.getChannelIndex() != a))
                      cout << "ERROR in col indices: " << endl;
                  } else if ((row == incol) && (col == inrow)) {
                    if ((K.getRowL() != L) || (K.getRowStimestwo() != S2) ||
                        (K.getRowChannelIndex() != a) ||
                        (K.getColumnL() != Lp) ||
                        (K.getColumnStimestwo() != S2p) ||
                        (K.getColumnChannelIndex() != ap))
                      cout << "ERROR in reversed indices" << endl;
                    if ((row.getL() != L) || (row.getStimestwo() != S2) ||
                        (row.getChannelIndex() != a))
                      cout << "ERROR in reversed row indices: " << endl;
                    if ((col.getL() != Lp) || (col.getStimestwo() != S2p) ||
                        (col.getChannelIndex() != ap))
                      cout << "ERROR in reversed col indices: " << endl;
                  } else
                    cout << "ERROR in row-col" << endl;
                  cout << K.str() << endl;
                  count1++;
                } catch (std::exception& xp) {
                  count2++;
                }
              }

  cout << "count1 = " << count1 << endl;
  cout << "count2 = " << count2 << endl;

  XMLHandler xmla(xmlin, "TestA");
  KIndex KK(xmla);
  cout << "KK = " << KK.str() << endl;

  XMLHandler xmlb(xmlin, "TestB");
  cout << xmlb.output() << endl;

  KElementInfo Kinfo(xmlb);
  cout << Kinfo.str() << endl;

  cout << endl << endl << "Testing KFitParamInfo" << endl << endl;
  string oname;

  KFitParamInfo KPtest;

  KFitParamInfo KP1;
  cout << "KP1 = " << KP1.output() << endl;
  oname = KP1.getMCObsName();
  cout << "KP1 string = " << oname << endl;
  KPtest.setFromMCObsName(oname);
  if (KPtest != KP1)
    cout << "ERROR" << endl;
  else
    cout << "Good read" << endl;

  KFitParamInfo KP2(Kinfo, 3);
  cout << "KP2 = " << KP2.output() << endl;
  oname = KP2.getMCObsName();
  cout << "KP2 string = " << oname << endl;
  KPtest.setFromMCObsName(oname);
  if (KPtest != KP2)
    cout << "ERROR" << endl;
  else
    cout << "Good read" << endl;

  KFitParamInfo KP3(2, 5);
  cout << "KP3 = " << KP3.output() << endl;
  oname = KP3.getMCObsName();
  cout << "KP3 string = " << oname << endl;
  KPtest.setFromMCObsName(oname);
  if (KPtest != KP3)
    cout << "ERROR" << endl;
  else
    cout << "Good read" << endl;

  KFitParamInfo KP4(KK, 1, 3);
  cout << "KP4 = " << KP4.output() << endl;
  oname = KP4.getMCObsName();
  cout << "KP4 string = " << oname << endl;
  KPtest.setFromMCObsName(oname);
  if (KPtest != KP4)
    cout << "ERROR" << endl;
  else
    cout << "Good read" << endl;

  XMLHandler xmlc(xmlin, "TestC");
  cout << xmlc.output() << endl;
  KFitParamInfo KP5(xmlc);
  cout << "KP5 = " << KP5.output() << endl;
  oname = KP5.getMCObsName();
  cout << "KP5 string = " << oname << endl;
  KPtest.setFromMCObsName(oname);
  if (KPtest != KP5)
    cout << "ERROR" << endl;
  else
    cout << "Good read" << endl;
}

void testKFitParamInfo(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestKFitParamInfo") == 0)
    return;

  cout << endl
       << endl
       << "***************************************************" << endl
       << endl;
  cout << "Testing KFitParamInfo" << endl;
  do_Kfitparamtests(xml_in);
}
