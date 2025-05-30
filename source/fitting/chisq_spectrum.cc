// #include "chisq_spectrum.h"
// using namespace std;

// SpectrumFit::~SpectrumFit() {
//   clear();
// }

// void SpectrumFit::clear() {
//   for (uint k = 0; k < BQ.size(); ++k)
//     delete BQ[k];
//   delete Kmat;
//   delete Kinv;
//   BQ.clear();
// }

// void SpectrumFit::guessInitialFitParamValues(
//     vector<double>& fitparams) const {
//   if (Kmat != 0)
//     fitparams = Kmat->getParameterValues();
//   else
//     fitparams = Kinv->getParameterValues();
// }