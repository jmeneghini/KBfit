#include "chisq_detres.h"
#include "xml_handler.h"
using namespace std;

void testChisqDetRes(XMLHandler& xml_in, int taskcount) {
  if (xml_tag_count(xml_in, "TestDeterminantResidualFit") == 0)
    return;

  cout << endl << "Starting TestDeterminantResidualFit" << endl;

  XMLHandler xmlr(xml_in, "DeterminantResidualFit");
  MCSamplingInfo sampinfo; // default jackknife
  KBObsHandler KBOH(sampinfo);
  XMLHandler xmlout;
  string outfileStub("fit_results");
  DeterminantResidualFit DRF(xmlr, &KBOH, xmlout, outfileStub);
  cout << endl << endl << DRF.output() << endl << endl;

  vector<double> fitparams;
  DRF.guessInitialFitParamValues(fitparams);

  LowerTriangularMatrix<double> inv_cov_cholesky;
  for (uint sindex = 0; sindex <= DRF.getNumberOfResamplings(); ++sindex) {
    cout << "sampling index = " << sindex << endl;
    DRF.setResamplingIndex(sindex);
    double chisq = DRF.evalChiSquare(fitparams);
    cout << "chi square = " << chisq << endl;
  }
}
