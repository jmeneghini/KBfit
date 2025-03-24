#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

std::string make_stringtuple(const std::vector<int>& ivals) {
  if (ivals.empty())
    return std::string("()");
  std::stringstream s;
  s << "(" << ivals[0];
  for (unsigned int k = 1; k < ivals.size(); ++k)
    s << "," << ivals[k];
  s << ")";
  return s.str();
}

std::vector<int> read_stringtuple(const std::string& istrlet) {
  std::vector<int> ivals;
  size_t pos1 = istrlet.find('(');
  size_t pos2 = istrlet.find(')', pos1 + 1);
  if ((pos1 == string::npos) || (pos2 == string::npos)) {
    return ivals;
  }
  for (size_t p = pos1 + 1; p < pos2;
       ++p) { // check for white space (not allowed)
    if (std::isspace(istrlet[p]))
      return ivals;
  }
  int k;
  size_t pos = istrlet.find_first_of(",)", pos1 + 1);
  while ((pos != string::npos) && (pos <= pos2)) {
    int count =
        sscanf(istrlet.substr(pos1 + 1, pos - pos1 - 1).c_str(), "%d", &k);
    if (count != 1) {
      ivals.clear();
      return ivals;
    }
    ivals.push_back(k);
    pos1 = pos;
    pos = istrlet.find_first_of(",)", pos1 + 1);
  }
  return ivals;
}

std::string make_stringtuple(int i1) {
  std::vector<int> ivals(1);
  ivals[0] = i1;
  return make_stringtuple(ivals);
}

std::string make_stringtuple(int i1, int i2) {
  std::vector<int> ivals(2);
  ivals[0] = i1;
  ivals[1] = i2;
  return make_stringtuple(ivals);
}

std::string make_stringtuple(int i1, int i2, int i3) {
  std::vector<int> ivals(3);
  ivals[0] = i1;
  ivals[1] = i2;
  ivals[2] = i3;
  return make_stringtuple(ivals);
}

void read_stringtuple(const std::string& str, int& i1) {
  std::vector<int> ivals = read_stringtuple(str);
  if (ivals.size() != 1)
    throw(std::invalid_argument("Bad stringtuple read"));
  i1 = ivals[0];
}

void read_stringtuple(const std::string& str, int& i1, int& i2) {
  std::vector<int> ivals = read_stringtuple(str);
  if (ivals.size() != 2)
    throw(std::invalid_argument("Bad stringtuple read"));
  i1 = ivals[0];
  i2 = ivals[1];
}

void read_stringtuple(const std::string& str, int& i1, int& i2, int& i3) {
  std::vector<int> ivals = read_stringtuple(str);
  if (ivals.size() != 3)
    throw(std::invalid_argument("Bad stringtuple read"));
  i1 = ivals[0];
  i2 = ivals[1];
  i3 = ivals[2];
}

std::vector<std::string> extract_stringtuples(const std::string& str) {
  std::vector<std::string> tuples;
  size_t pos = str.find('(');
  while (pos != string::npos) {
    cout << "pos = " << pos << endl;
    size_t pos1 = str.find('(', pos + 1);
    cout << "pos1 = " << pos1 << endl;
    size_t pos2 = str.find(')', pos + 1);
    cout << "pos2 = " << pos2 << endl;
    if ((pos2 != string::npos) && (pos1 > pos2))
      tuples.push_back(str.substr(pos, pos2 - pos + 1));
    pos = str.find('(', pos2);
    cout << "pos=" << pos << endl;
  }
  return tuples;
}

int main() {

  cout << make_stringtuple(7) << endl;
  cout << make_stringtuple(3, 4) << endl;
  cout << make_stringtuple(0, 2) << endl;
  cout << make_stringtuple(0, 0) << endl;
  cout << make_stringtuple(-1, 4) << endl;
  cout << make_stringtuple(3, -4, 6) << endl;

  string tester("   (19,0,8) ");
  vector<int> ivals = read_stringtuple(tester);
  cout << "size of ivals = " << ivals.size() << endl;
  for (unsigned int k = 0; k < ivals.size(); ++k)
    cout << "ivals[" << k << "] = " << ivals[k] << endl;

  int i1, i2, i3;
  read_stringtuple(tester, i1, i2, i3);
  cout << "i1 = " << i1 << endl;
  cout << "i2 = " << i2 << endl;
  cout << "i3 = " << i3 << endl;

  std::vector<std::string> tp =
      extract_stringtuples("  adfdaf a (4,5)  (() (6,5)");
  for (unsigned int k = 0; k < tp.size(); ++k)
    cout << tp[k] << endl;

  return 0;
}
