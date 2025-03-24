#include "ni_reg.h"
#include "xml_handler.h"
using namespace std;

void testNonIntRegulator(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestNonIntRegulator") == 0)
    return;

  cout << "TestNonIntRegulator" << endl;

  vector<double> s(3, 0.0);
  double gam = 1.0;
  double usq_min = 0.0;
  double usq_max = 11.6;

  NonIntRegulator NIR(s, gam, usq_min, usq_max);
  double usq = usq_min;
  while (usq <= usq_max) {
    cout << " " << usq << "   " << NIR.getValue(usq) << endl;
    usq += 0.01;
  }
}

// ******************************************************************************
