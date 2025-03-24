#ifndef KBFIT_TESTS_H
#define KBFIT_TESTS_H

#include "xml_handler.h"

// ************************************************

void testMatrix(XMLHandler& xmlr);
// void testBootstrapper(XMLHandler& xmlr);
// void testCorrDataHandler(XMLHandler& xmlr);
// void testCorrelatorInfo(XMLHandler& xmlr);
// void testCorrMatEstimates(XMLHandler& xmlr);
// void testVEVDataHandler(XMLHandler& xmlr);
// void testTaskHandler(XMLHandler& xmlr);
// void testXMLHandler(XMLHandler& xmlr);
// void testArgsHandler(XMLHandler& xmlr);
// void testOperatorInfo(XMLHandler& xmlr);
void testMCObsInfo(XMLHandler& xmlr);
void testKBObsInfo(XMLHandler& xmlr);
void testMCObsGetHandler(XMLHandler& xml_in);
void testMCObsGetHandlerFake(XMLHandler& xml_in);
// void testMCObsHandler(XMLHandler& xml_in);
// void testMCObsHandler2(XMLHandler& xml_in);
// void testMCObsHandlerIO(XMLHandler& xml_in);
// void testGracePlot(XMLHandler& xml_in);
// void testmulticompare(XMLHandler& xml_in);
// void testCorrelatorMatrixInfo(XMLHandler& xml_in);
// void testChiSquare(XMLHandler& xml_in, int taskcount);
// void testGamma(XMLHandler& xml_in);
// void testEffEnergy(XMLHandler& xml_in);
void testChisqDetRes(XMLHandler& xml_in, int taskcount);
// void testRotateCorrelator(XMLHandler& xml_in, int taskcount);
// void testPivotCorrelator(XMLHandler& xml_in, int taskcount);
// void testPivotCorrelator0(XMLHandler& xml_in, int taskcount);
void testMCBinsInfo(XMLHandler& xml_in);
void testMCSamplingInfo(XMLHandler& xml_in);
void testMinimizer(XMLHandler& xml_in);
void testKFitParamInfo(XMLHandler& xml_in);
// void testDataGetPut(XMLHandler& xml_in);
// void testBinGetPut(XMLHandler& xml_in);
// void testSamplingGetPut(XMLHandler& xml_in);
// void testReorder(XMLHandler& xml_in);
// void testNonIntRegulator(XMLHandler& xml_in);
void testMatrixInverse(XMLHandler& xml_in);

// ***********************************************
#endif
