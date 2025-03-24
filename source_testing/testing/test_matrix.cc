#include "matrix.h"
#include "xml_handler.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <list>
#include <set>
using namespace std;

//   outvec = inmat * invec

void multiply(CVector& outvec, const ComplexHermitianMatrix& inmat,
              const CVector& invec) {
  uint nn = inmat.size();
  if (nn != invec.size())
    throw(std::invalid_argument("bad matrix-vector multiplying sizes"));
  outvec.resize(nn);
  for (uint row = 0; row < nn; row++) {
    complex<double> z(0.0, 0.0);
    for (uint k = 0; k < nn; k++)
      z += inmat(row, k) * invec[k];
    outvec[row] = z;
  }
}

void multiply(RVector& outvec, const RealSymmetricMatrix& inmat,
              const RVector& invec) {
  uint nn = inmat.size();
  if (nn != invec.size())
    throw(std::invalid_argument("bad matrix-vector multiplying sizes"));
  outvec.resize(nn);
  for (uint row = 0; row < nn; row++) {
    double z = 0.0;
    for (uint k = 0; k < nn; k++)
      z += inmat(row, k) * invec[k];
    outvec[row] = z;
  }
}

//   inner product of vectors

complex<double> dotProduct(const CVector& lvec, const CVector& rvec) {
  uint nv = rvec.size();
  if (nv != lvec.size())
    throw(std::invalid_argument("Mismatch in vector sizes in dotProduct"));
  complex<double> res(0.0, 0.0);
  for (uint k = 0; k < nv; k++)
    res += conjugate(lvec[k]) * rvec[k];
  return res;
}

double dotProduct(const RVector& lvec, const RVector& rvec) {
  uint nv = rvec.size();
  if (nv != lvec.size())
    throw(std::invalid_argument("Mismatch in vector sizes in dotProduct"));
  double res = 0.0;
  for (uint k = 0; k < nv; k++)
    res += lvec[k] * rvec[k];
  return res;
}

// *******************************************************************************
// * *
// *    The classes "MultiIntLooper" and "MultiIntLooperWR" (with restrictions)
// *
// *    are defined in this file.  These are used in the CorrelatorHandler *
// *    compute routines to loop over dilution indices and noise choices. *
// * *
// *    The following utility routines are also defined here: *
// *          "make_list",  "make_vector",   "cmplx_same" *
// * *
// *******************************************************************************

//   Simple routine to help make a list of objects

template <typename T> std::list<T> make_list() {
  std::list<T> result;
  return result;
}

template <typename T> std::list<T> make_list(const T& v1) {
  std::list<T> result;
  result.push_back(v1);
  return result;
}

template <typename T> std::list<T> make_list(const T& v1, const T& v2) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  result.push_back(v4);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  result.push_back(v4);
  result.push_back(v5);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  result.push_back(v4);
  result.push_back(v5);
  result.push_back(v6);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  result.push_back(v4);
  result.push_back(v5);
  result.push_back(v6);
  result.push_back(v7);
  return result;
}

template <typename T>
std::list<T> make_list(const T& v1, const T& v2, const T& v3, const T& v4,
                       const T& v5, const T& v6, const T& v7, const T& v8) {
  std::list<T> result;
  result.push_back(v1);
  result.push_back(v2);
  result.push_back(v3);
  result.push_back(v4);
  result.push_back(v5);
  result.push_back(v6);
  result.push_back(v7);
  result.push_back(v8);
  return result;
}

//   Simple routine to help make a vector of objects

template <typename T> std::vector<T> make_vector() {
  std::vector<T> result;
  return result;
}

template <typename T> std::vector<T> make_vector(const T& v1) {
  std::vector<T> result(1);
  result[0] = v1;
  return result;
}

template <typename T> std::vector<T> make_vector(const T& v1, const T& v2) {
  std::vector<T> result(2);
  result[0] = v1;
  result[1] = v2;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3) {
  std::vector<T> result(3);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4) {
  std::vector<T> result(4);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  result[3] = v4;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5) {
  std::vector<T> result(5);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  result[3] = v4;
  result[4] = v5;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6) {
  std::vector<T> result(6);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  result[3] = v4;
  result[4] = v5;
  result[5] = v6;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7) {
  std::vector<T> result(7);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  result[3] = v4;
  result[4] = v5;
  result[5] = v6;
  result[6] = v7;
  return result;
}

template <typename T>
std::vector<T> make_vector(const T& v1, const T& v2, const T& v3, const T& v4,
                           const T& v5, const T& v6, const T& v7, const T& v8) {
  std::vector<T> result(8);
  result[0] = v1;
  result[1] = v2;
  result[2] = v3;
  result[3] = v4;
  result[4] = v5;
  result[5] = v6;
  result[6] = v7;
  result[7] = v8;
  return result;
}

// ********************************************************************************
// * *
// *    Multi-dimensional looper over unsigned integers.  For dimension "k", *
// *    its index varies from 0..upper[k]-1.  The number of dimensions *
// *    and the upper limits must be specified in the constructor.  Each *
// *    dimension must have value 1 or larger (zero values are not allowed). *
// * *
// *    Usage:  MultiIntLooper v(3,4);   // row value 0,1,2   column = 0,1,2,3 *
// *            for (v.start();v.notdone();++v){....} *
// *     loop order = v(0,0), v(1,0), v(2,0), v(0,1), v(1,1), v(2,1), ... *
// * *
// ********************************************************************************

class MultiIntLooper {

  unsigned int m_dim;
  std::vector<unsigned int> m_upper;
  std::vector<unsigned int> m_current;
  bool m_notdone;

public:
  MultiIntLooper(const std::vector<unsigned int> upper_limits)
      : m_dim(upper_limits.size()), m_upper(upper_limits),
        m_current(upper_limits.size(), 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(unsigned int upper1)
      : m_dim(1), m_upper(make_vector(upper1)), m_current(1, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(int upper1)
      : m_dim(1), m_upper(make_vector(unsigner(upper1))), m_current(1, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(unsigned int upper1, unsigned int upper2)
      : m_dim(2), m_upper(make_vector(upper1, upper2)), m_current(2, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(int upper1, int upper2)
      : m_dim(2), m_upper(make_vector(unsigner(upper1), unsigner(upper2))),
        m_current(2, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(unsigned int upper1, unsigned int upper2, unsigned int upper3)
      : m_dim(3), m_upper(make_vector(upper1, upper2, upper3)), m_current(3, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(int upper1, int upper2, int upper3)
      : m_dim(3), m_upper(make_vector(unsigner(upper1), unsigner(upper2),
                                      unsigner(upper3))),
        m_current(3, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(unsigned int upper1, unsigned int upper2, unsigned int upper3,
                 unsigned int upper4)
      : m_dim(4), m_upper(make_vector(upper1, upper2, upper3, upper4)),
        m_current(4, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooper(int upper1, int upper2, int upper3, int upper4)
      : m_dim(4), m_upper(make_vector(unsigner(upper1), unsigner(upper2),
                                      unsigner(upper3), unsigner(upper4))),
        m_current(4, 0), m_notdone(true) {
    checkUpper();
  }

  ~MultiIntLooper() {}

  const MultiIntLooper& start() {
    for (unsigned int k = 0; k < m_dim; ++k)
      m_current[k] = 0;
    m_notdone = true;
    return *this;
  }

  bool notdone() const { return m_notdone; }

  const MultiIntLooper& operator++() // prefix
  {
    return increment();
  }

  unsigned int numDimensions() const { return m_dim; }

  unsigned int operator[](unsigned int k) const { return m_current[k]; }

  const std::vector<unsigned int>& get_current() const { return m_current; }

private:
  const MultiIntLooper& increment() {
    if (!m_notdone)
      return *this;
    ++m_current[0];
    unsigned int j = 0;
    while (m_current[j] == m_upper[j]) {
      m_current[j] = 0;
      if (j == (m_dim - 1)) {
        m_notdone = false;
        break;
      }
      ++m_current[++j];
    }
    return *this;
  }

  void checkUpper() {
    if (m_upper.size() == 0)
      throw(std::invalid_argument("Bad MultiIntLooper: no dimensions"));
    for (unsigned int k = 0; k < m_dim; ++k)
      if (m_upper[k] < 1)
        throw(std::invalid_argument("Bad MultiIntLooper upper limits"));
  }

  unsigned int unsigner(int k) { return (k >= 0) ? (unsigned int)(k) : 0; }
};

// ********************************************************************************
// * *
// *    Multi-dimensional looper over unsigned integers.  For dimension "k", *
// *    its index varies from 0..upper[k]-1.  The number of dimensions *
// *    and the upper limits must be specified in the constructor.  Each *
// *    dimension must have value 1 or larger (zero values are not allowed). *
// *    The "WR" here means "with restrictions".  This looper will avoid any *
// *    instances in which any two indices have the same values, except for *
// *    those pairs explicitly allowed by an "allowEqual" call. *
// * *
// *    Usage:  MultiIntLooperWR v(3,4);   // row value 0,1,2   column = 0,1,2,3
// *
// *            v.allowEqual(0,1);          // index 0 and 1 can be equal *
// *            for (v.start();v.notdone();++v){....} *
// *     loop order = v(0,0), v(1,0), v(2,0), v(0,1), v(1,1), v(2,1), ... *
// * *
// ********************************************************************************

class MultiIntLooperWR {

  unsigned int m_dim;
  std::vector<unsigned int> m_upper;
  std::vector<unsigned int> m_current;
  bool m_notdone;
  std::set<std::pair<unsigned int, unsigned int>> m_allowEqual;

public:
  MultiIntLooperWR(const std::vector<unsigned int> upper_limits)
      : m_dim(upper_limits.size()), m_upper(upper_limits),
        m_current(upper_limits.size(), 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(unsigned int upper1)
      : m_dim(1), m_upper(make_vector(upper1)), m_current(1, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(int upper1)
      : m_dim(1), m_upper(make_vector(unsigner(upper1))), m_current(1, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(unsigned int upper1, unsigned int upper2)
      : m_dim(2), m_upper(make_vector(upper1, upper2)), m_current(2, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(int upper1, int upper2)
      : m_dim(2), m_upper(make_vector(unsigner(upper1), unsigner(upper2))),
        m_current(2, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(unsigned int upper1, unsigned int upper2,
                   unsigned int upper3)
      : m_dim(3), m_upper(make_vector(upper1, upper2, upper3)), m_current(3, 0),
        m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(int upper1, int upper2, int upper3)
      : m_dim(3), m_upper(make_vector(unsigner(upper1), unsigner(upper2),
                                      unsigner(upper3))),
        m_current(3, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(unsigned int upper1, unsigned int upper2,
                   unsigned int upper3, unsigned int upper4)
      : m_dim(4), m_upper(make_vector(upper1, upper2, upper3, upper4)),
        m_current(4, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR(int upper1, int upper2, int upper3, int upper4)
      : m_dim(4), m_upper(make_vector(unsigner(upper1), unsigner(upper2),
                                      unsigner(upper3), unsigner(upper4))),
        m_current(4, 0), m_notdone(true) {
    checkUpper();
  }

  MultiIntLooperWR& allowEqual(int index1, int index2) {
    return allowEqual((unsigned int)(index1), (unsigned int)(index2));
  }

  MultiIntLooperWR& allowEqual(unsigned int index1, int unsigned index2) {
    if ((index1 < index2) && (index2 < m_dim))
      m_allowEqual.insert(std::make_pair(index1, index2));
    else if ((index2 > index1) && (index1 < m_dim))
      m_allowEqual.insert(std::make_pair(index2, index1));
    return *this;
  }

  ~MultiIntLooperWR() {}

  const MultiIntLooperWR& start() {
    for (unsigned int k = 0; k < m_dim; ++k)
      m_current[k] = 0;
    m_notdone = true;
    if (m_dim < 2)
      return *this;
    if (!allowed())
      operator++();
    return *this;
  }

  bool notdone() const { return m_notdone; }

  const MultiIntLooperWR& operator++() // prefix
  {
    increment();
    while ((!allowed()) && (m_notdone))
      increment();
    return *this;
  }

  unsigned int numDimensions() const { return m_dim; }

  unsigned int operator[](unsigned int k) const { return m_current[k]; }

  const std::vector<unsigned int>& get_current() const { return m_current; }

private:
  const MultiIntLooperWR& increment() {
    if (!m_notdone)
      return *this;
    ++m_current[0];
    unsigned int j = 0;
    while (m_current[j] == m_upper[j]) {
      m_current[j] = 0;
      if (j == (m_dim - 1)) {
        m_notdone = false;
        break;
      }
      ++m_current[++j];
    }
    return *this;
  }

  bool allowed() {
    for (unsigned int j = 0; j < (m_dim - 1); ++j)
      for (unsigned int k = j + 1; k < m_dim; ++k)
        if ((m_current[j] == m_current[k]) &&
            (m_allowEqual.count(std::make_pair(j, k)) == 0))
          return false;
    return true;
  }

  void checkUpper() {
    if (m_upper.size() == 0)
      throw(std::invalid_argument("Bad MultiIntLooperWR: no dimensions"));
    for (unsigned int k = 0; k < m_dim; ++k)
      if (m_upper[k] < 1)
        throw(std::invalid_argument("Bad MultiIntLooperWR upper limits"));
  }

  unsigned int unsigner(int k) { return (k >= 0) ? (unsigned int)(k) : 0; }
};

// ****************************************************************************

inline bool equal(const float& x1, const float& x2) {
  float abstol = 1e-8, reltol = 1e-5;
  return (abs(x1 - x2) <= (abstol + reltol * abs(x2)));
}

inline bool equal(const std::complex<float>& z1,
                  const std::complex<float>& z2) {
  return (equal(real(z1), real(z2)) && equal(imag(z1), imag(z2)));
}

inline bool equal(const double& x1, const double& x2) {
  double abstol = 1e-15, reltol = 1e-12;
  return (abs(x1 - x2) <= (abstol + reltol * abs(x2)));
}

inline bool equal(const std::complex<double>& z1,
                  const std::complex<double>& z2) {
  return (equal(real(z1), real(z2)) && equal(imag(z1), imag(z2)));
}

inline bool equal(int a, int b) { return (a == b); }

inline bool equal(unsigned int a, unsigned int b) { return (a == b); }

unsigned int inc(unsigned int& k, unsigned int maximum) {
  unsigned int v = k;
  k++;
  if (k == maximum)
    k = 0;
  return v;
}

// ****************************************************************************

void dotestVectorSimple() {
  cout << endl << endl << "Simple Test Results" << endl << endl;

  RVector A0;
  cout << "A0 size should be 0:  result = " << A0.size() << endl;
  try {
    cout << "Trying A0[5]: should throw exception" << endl;
    cout << A0[5] << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  int ival = 4;
  RVector A1(ival);
  cout << "A1 size should be 4:  result = " << A1.size() << endl;
  uint uival = 6;
  RVector A1b(uival);
  cout << "A1b size should be 6:  result = " << A1b.size() << endl;

  ival = 3;
  RVector A2(ival, 4.6);
  cout << "A2 size should be 3:  result = " << A2.size() << endl;
  for (int k = 0; k < int(A2.size()); ++k) {
    cout << "Should be initialized to 4.6:  result[" << k << "] = " << A2[k]
         << endl;
    cout << "Should be initialized to 4.6:  result(" << k << ") = " << A2(k)
         << endl;
  }
  for (uint k = 0; k < A2.size(); ++k) {
    cout << "Should be initialized to 4.6:  result[" << k << "] = " << A2[k]
         << endl;
    cout << "Should be initialized to 4.6:  result(" << k << ") = " << A2(k)
         << endl;
  }
  A2[2] = 1.5;
  cout << "Should be assigned 1.5:  result[" << 2 << "] = " << A2[2] << endl;
  A2(1) = 11.5;
  cout << "Should be assigned 11.5:  result[" << 1 << "] = " << A2(1) << endl;
  A2.put(0, 54.3);
  cout << "Should be 54.3 put:  result[" << 0 << "] = " << A2.get(0) << endl;

  uival = 5;
  RVector A2b(uival, 14.6);
  cout << "A2b size should be 5:  result = " << A2b.size() << endl;
  for (int k = 0; k < int(A2b.size()); ++k) {
    cout << "Should be initialized to 14.6:  result[" << k << "] = " << A2b[k]
         << endl;
    cout << "Should be initialized to 14.6:  result(" << k << ") = " << A2b(k)
         << endl;
  }
  for (uint k = 0; k < A2b.size(); ++k) {
    cout << "Should be initialized to 14.6:  result[" << k << "] = " << A2b[k]
         << endl;
    cout << "Should be initialized to 14.6:  result(" << k << ") = " << A2b(k)
         << endl;
  }
  A2[2] = 1.5;
  cout << "Should be assigned 1.5:  result[" << 2 << "] = " << A2[2] << endl;
  A2(1) = 11.5;
  cout << "Should be assigned 11.5:  result[" << 1 << "] = " << A2(1) << endl;
  A2.put(0, 54.3);
  cout << "Should be 54.3 put:  result[" << 0 << "] = " << A2.get(0) << endl;

  ival = 8;
  A2.resize(ival);
  cout << "A2 should be resized to 8: size = " << A2.size() << endl;
  uival = 7;
  A2.put(uival, 36.7);
  cout << "Should be 36.7 put:  result[" << 7 << "] = " << A2.get(7) << endl;
  try {
    cout << "Trying A2[12]: should throw exception" << endl;
    cout << A2[12] << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  A2.resize();
  cout << "A2 should be resized to 0: size = " << A2.size() << endl;
  try {
    cout << "Trying to resize to -3: should throw exception" << endl;
    A2.resize(-3);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  A1.clear();
  cout << "A1 cleared: size should be 0:  result = " << A1.size() << endl;

  uival = 5;
  A0.resize(uival);
  A0[0] = 2.3;
  A0[1] = 33.0;
  A0[2] = -5.6;
  A0[3] = 88.8;
  A0[4] = 10.0;
  RVector A3(A0);
  cout << "size of A3 should be 5: size = " << A3.size() << endl;
  cout << "A3[0] should be 2.3:  result = " << A3[0] << endl;
  cout << "A3[1] should be 33.0:  result = " << A3[1] << endl;
  cout << "A3[2] should be -5.6:  result = " << A3[2] << endl;
  cout << "A3[3] should be 88.8:  result = " << A3[3] << endl;
  cout << "A3[4] should be 10:  result = " << A3[4] << endl;

  vector<double> stdvec(4);
  stdvec[0] = 2.27;
  stdvec[1] = 3.3;
  stdvec[2] = -5.6;
  stdvec[3] = 7.4;
  RVector A4(stdvec);
  cout << "size of A4 should be 4: size = " << A4.size() << endl;
  cout << "A4[0] should be 2.27:  result = " << A4[0] << endl;
  cout << "A4[1] should be 3.3:  result = " << A4[1] << endl;
  cout << "A4[2] should be -5.6:  result = " << A4[2] << endl;
  cout << "A4[3] should be 7.4:  result = " << A4[3] << endl;

  A1.resize(5);
  A1 = 7.7;
  cout << "A1 should now be size 5: size = " << A1.size() << endl;
  for (int k = 0; k < int(A1.size()); ++k) {
    cout << "Should be set to 7.7:  result[" << k << "] = " << A1[k] << endl;
  }
  A1 += 3.3;
  for (int k = 0; k < int(A1.size()); ++k) {
    cout << "Should be += to 11.0:  result[" << k << "] = " << A1[k] << endl;
  }
  A1 -= 6.0;
  for (int k = 0; k < int(A1.size()); ++k) {
    cout << "Should be -= to 5.0:  result[" << k << "] = " << A1[k] << endl;
  }
  A1 *= 3.0;
  for (int k = 0; k < int(A1.size()); ++k) {
    cout << "Should be *= to 15.0:  result[" << k << "] = " << A1[k] << endl;
  }
  A1 /= 2.0;
  for (int k = 0; k < int(A1.size()); ++k) {
    cout << "Should be /= to 7.5:  result[" << k << "] = " << A1[k] << endl;
  }

  A2 = A4;
  cout << "A2=A4:  size of A2 should be 4: size = " << A2.size() << endl;
  cout << "A2=A4:  A2[0] should be 2.27:  result = " << A2[0] << endl;
  cout << "A2=A4:  A2[1] should be 3.3:  result = " << A2[1] << endl;
  cout << "A2=A4:  A2[2] should be -5.6:  result = " << A2[2] << endl;
  cout << "A2=A4:  A2[3] should be 7.4:  result = " << A2[3] << endl;

  try {
    cout << "Trying to add mismatched size: should throw exception" << endl;
    A1 += A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to subtract mismatched size: should throw exception"
         << endl;
    A1 -= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to multiply mismatched size: should throw exception"
         << endl;
    A1 *= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to divide mismatched size: should throw exception" << endl;
    A1 /= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }

  A1 = A2;
  A2 += A1;
  cout << "A2+=A1:  size of A2 should be 4: size = " << A2.size() << endl;
  cout << "A2+=A1:  A2[0] should be 4.54.:  result = " << A2[0] << endl;
  cout << "A2+=A1:  A2[1] should be 6.6:  result = " << A2[1] << endl;
  cout << "A2+=A1:  A2[2] should be -11.2:  result = " << A2[2] << endl;
  cout << "A2+=A1:  A2[3] should be 14.8:  result = " << A2[3] << endl;

  A1[0] = 1.1;
  A1[1] = 2.2;
  A1[2] = 3.3;
  A1[3] = 4.4;
  A2 -= A1;
  cout << "A2-=A1:  size of A2 should be 4: size = " << A2.size() << endl;
  cout << "A2-=A1:  A2[0] should be 3.44.:  result = " << A2[0] << endl;
  cout << "A2-=A1:  A2[1] should be 4.4:  result = " << A2[1] << endl;
  cout << "A2-=A1:  A2[2] should be -14.5:  result = " << A2[2] << endl;
  cout << "A2-=A1:  A2[3] should be 10.4:  result = " << A2[3] << endl;

  A2 *= A2;
  cout << "A2*=A2:  size of A2 should be 4: size = " << A2.size() << endl;
  cout << "A2*=A2:  A2[0] should be 11.8336:  result = " << A2[0] << endl;
  cout << "A2*=A2:  A2[1] should be 19.36:  result = " << A2[1] << endl;
  cout << "A2*=A2:  A2[2] should be 210.25:  result = " << A2[2] << endl;
  cout << "A2*=A2:  A2[3] should be 108.16:  result = " << A2[3] << endl;

  A2 /= A1;
  cout << "A2/=A1:  size of A2 should be 4: size = " << A2.size() << endl;
  cout << "A2/=A1:  A2[0] should be 10.757818182:  result = " << A2[0] << endl;
  cout << "A2/=A1:  A2[1] should be 8.8:  result = " << A2[1] << endl;
  cout << "A2/=A1:  A2[2] should be 63.712121212:  result = " << A2[2] << endl;
  cout << "A2/=A1:  A2[3] should be 24.581818182:  result = " << A2[3] << endl;

  vector<double> res = A2.c_vector();
  cout << "res=c_vector:  size of res should be 4: size = " << res.size()
       << endl;
  cout << "res=c_vector:  res[0] should be 10.757818182:  result = " << res[0]
       << endl;
  cout << "res=c_vector:  res[1] should be 8.8:  result = " << res[1] << endl;
  cout << "res=c_vector:  res[2] should be 63.712121212:  result = " << res[2]
       << endl;
  cout << "res=c_vector:  res[3] should be 24.581818182:  result = " << res[3]
       << endl;

  cout << "simple tests done" << endl << endl;
}

template <typename T>
bool dotestVector(const vector<T>& datacmp, const vector<unsigned int>& idata) {
  unsigned int successcount = 0;
  unsigned int failcount = 0;
  unsigned int datacount = 0;
  unsigned int icount = 0;
  unsigned int dmax = datacmp.size();
  unsigned int imax = idata.size();

  unsigned int ndim = 1;
  Vector<T>* Aptr = 0;
  vector<unsigned int> csizes;
  Aptr = new Vector<T>(2);
  csizes.resize(ndim);
  csizes[0] = 2;

  unsigned int sizecheck = Aptr->size();
  unsigned int ss = 1;
  for (uint k = 0; k < ndim; ++k)
    ss *= csizes[k];
  if (ss == sizecheck)
    successcount++;
  else {
    failcount++;
    cout << "size test failed" << endl;
  }
  delete Aptr;
  Aptr = 0;

  csizes.resize(ndim);
  for (uint i = 0; i < ndim; ++i)
    csizes[i] = idata[inc(icount, imax)];
  Vector<T> A0(csizes[0]);

  for (unsigned int jj = 0; jj < 3; jj++) {

    if (jj != 0) {
      for (uint i = 0; i < ndim; ++i)
        csizes[i] = idata[inc(icount, imax)];
      A0.resize(csizes[0]);
    }

    unsigned int sizecheck = A0.size();
    unsigned int ss = 1;
    for (uint k = 0; k < ndim; ++k)
      ss *= csizes[k];
    if (ss == sizecheck)
      successcount++;
    else {
      failcount++;
      cout << "size test failed" << endl;
    }

    datacount = 0;
    unsigned int dmax = datacmp.size();
    for (uint k0 = 0; k0 < csizes[0]; ++k0) {
      A0(k0) = datacmp[inc(datacount, dmax)];
      if (datacount >= dmax)
        datacount = 0;
    }

    datacount = 0;
    for (uint k0 = 0; k0 < csizes[0]; ++k0) {
      if (A0(k0) != datacmp[inc(datacount, dmax)])
        failcount++;
      else
        successcount++;
      if (datacount >= dmax)
        datacount = 0;
    }

    datacount = 0;
    dmax = datacmp.size();
    for (uint k0 = 0; k0 < csizes[0]; ++k0) {
      A0.put(k0, datacmp[inc(datacount, dmax)]);
      if (datacount >= dmax)
        datacount = 0;
    }

    datacount = 0;
    for (uint k0 = 0; k0 < csizes[0]; ++k0) {
      if (A0.get(k0) != datacmp[inc(datacount, dmax)])
        failcount++;
      else
        successcount++;
      if (datacount >= dmax)
        datacount = 0;
    }
  }

  vector<unsigned int> isizes(ndim);
  for (uint k = 0; k < ndim; ++k)
    isizes[k] = idata[inc(icount, imax)];
  T cf = datacmp[inc(datacount, dmax)];
  Vector<T> B0(isizes[0], cf);
  MultiIntLooper ind(isizes);
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  cf = datacmp[inc(datacount, dmax)];
  cf *= datacmp[inc(datacount, dmax)];
  cf -= datacmp[inc(datacount, dmax)];
  B0 = cf;
  cout << "cf = " << cf << endl;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  T cf2 = datacmp[inc(datacount, dmax)];
  cf += cf2;
  B0 += cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf *= cf2;
  B0 *= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf -= cf2;
  B0 -= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf /= cf2;
  B0 /= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()[0]), cf))
      successcount++;
    else
      failcount++;
  }

  Vector<T> W(B0);
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(W(ind.get_current()[0]), B0(ind.get_current()[0])))
      successcount++;
    else
      failcount++;
  }

  for (uint k = 0; k < ndim; ++k)
    isizes[k] = idata[inc(icount, imax)];
  B0.resize(isizes[0]);
  MultiIntLooper dil(isizes);
  datacount = 0;
  for (dil.start(); dil.notdone(); ++dil)
    B0(dil.get_current()[0]) = datacmp[inc(datacount, dmax)];
  datacount = dmax / 2;
  Vector<T> B1(isizes[0]);
  for (dil.start(); dil.notdone(); ++dil)
    B1(dil.get_current()[0]) = datacmp[inc(datacount, dmax)];

  Vector<T> B2(B0);
  B2 += B1;
  datacount = 0;
  unsigned int dcount2 = dmax / 2;
  for (dil.start(); dil.notdone(); ++dil) {
    if (equal(B2(dil.get_current()[0]),
              (datacmp[inc(datacount, dmax)] + datacmp[inc(dcount2, dmax)])))
      successcount++;
    else
      failcount++;
  }

  Vector<T> B3(B0);
  B3 *= B1;
  datacount = 0;
  dcount2 = dmax / 2;
  for (dil.start(); dil.notdone(); ++dil) {
    if (equal(B3(dil.get_current()[0]),
              (datacmp[inc(datacount, dmax)] * datacmp[inc(dcount2, dmax)])))
      successcount++;
    else
      failcount++;
  }

  cout << endl;
  cout << " Number of successful tests = " << successcount << endl;
  cout << "     Number of FAILED tests = " << failcount << endl;
  return (failcount == 0);
}

void testVector() {

  cout << "Starting Vector tests" << endl;

  dotestVectorSimple();

  bool success = true;
  srand(time(NULL));
  unsigned int nintegers = 2056;
  vector<unsigned int> idata(nintegers);
  for (unsigned int k = 0; k < nintegers; ++k)
    idata[k] = rand() % 8 + 1;

  unsigned int ndata = 5048;
  vector<int> intdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    intdata[k] = (rand() % 4096 + 1) - 2047;

  cout << endl << "Tests: integer type" << endl;
  success = success && dotestVector(intdata, idata);
  cout << endl;

  ndata = 5048;
  vector<unsigned int> uintdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    uintdata[k] = (rand() % 4096 + 1);

  cout << endl << "Tests: unsigned integer type" << endl;
  success = success && dotestVector(uintdata, idata);
  cout << endl;

  ndata = 5048;
  vector<float> fdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    fdata[k] =
        1.434 * ((rand() % 4096 + 1) - 2047) - 0.54511 * ((rand() % 128) - 63);

  cout << endl << "Tests: float type" << endl;
  success = success && dotestVector(fdata, idata);
  cout << endl;

  ndata = 5048;
  vector<double> ddata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    ddata[k] = -2.3418 * ((rand() % 4096 + 1) - 2047) -
               0.22485 * ((rand() % 128) - 63);

  cout << endl << "Tests: double type" << endl;
  success = success && dotestVector(ddata, idata);
  cout << endl;

  ndata = 5048;
  vector<complex<float>> zfdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    zfdata[k] = complex<float>(1.434 * ((rand() % 4096 + 1) - 2047),
                               -0.54511 * ((rand() % 128) - 63));

  cout << endl << "Tests: complex<float> type" << endl;
  success = success && dotestVector(zfdata, idata);
  cout << endl;

  ndata = 5048;
  vector<complex<double>> zddata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    zddata[k] = complex<double>(-2.3418 * ((rand() % 4096 + 1) - 2047),
                                0.22485 * ((rand() % 128) - 63));

  cout << endl << "Tests: complex<double> type" << endl;
  success = success && dotestVector(zddata, idata);
  cout << endl;

  vector<double> vlast(5);
  for (int i = 0; i < 5; i++)
    vlast[i] = i * 2.2;
  Vector<double> VL(vlast);
  bool vflag = true;
  for (int i = 0; i < 5; i++)
    if (std::abs(VL[i] - vlast[i]) > 1e-6)
      vflag = false;
  vector<double> vtemp(VL.c_vector());
  VL[0] = 12.0;
  for (int i = 0; i < 5; i++)
    if (std::abs(vtemp[i] - vlast[i]) > 1e-6)
      vflag = false;

  success = success && vflag;

  if (success)
    cout << "ALL VECTOR TESTS PASSED!!" << endl;
  else
    cout << "Some vector tests FAILED" << endl;
}

// *******************************************************************************

void dotestMatrixSimple() {
  cout << endl << endl << "Simple Test Results" << endl << endl;

  RMatrix A0;
  cout << "A0 size should be 0:  result = " << A0.size() << endl;
  try {
    cout << "Trying A0(2,3): should throw exception" << endl;
    cout << A0(2, 3) << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  int ival1 = 4, ival2 = 3;
  RMatrix A1(ival1, ival2);
  cout << "A1 size should be 12:  result = " << A1.size() << endl;
  ival1 = 0;
  cout << "A1.size(0) should be 4: result = " << A1.size(ival1) << endl;
  ival1 = 1;
  cout << "A1.size(1) should be 3: result = " << A1.size(ival1) << endl;
  uint uival1 = 0;
  cout << "A1.size(0) should be 4: result = " << A1.size(uival1) << endl;
  uival1 = 1;
  cout << "A1.size(1) should be 3: result = " << A1.size(uival1) << endl;
  try {
    cout << "A1.size(2) should throw exception " << A1.size(2) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  ival1 = 3;
  ival2 = 2;
  RMatrix A2(ival1, ival2, 4.6);
  cout << "A2 size should be 6:  result = " << A2.size() << endl;
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "Should be initialized to 4.6:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  try {
    cout << "A2(4,5) should throw exception" << endl;
    cout << A2(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2(0,5) should throw exception" << endl;
    cout << A2(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  uival1 = 4;
  uint uival2 = 3;
  RMatrix A2b(uival1, uival2, 24.6);
  cout << "A2b size should be 12:  result = " << A2b.size() << endl;
  for (uint j = 0; j < A2b.size(0); ++j)
    for (uint k = 0; k < A2b.size(1); ++k) {
      cout << "Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  try {
    cout << "A2b(4,5) should throw exception" << endl;
    cout << A2b(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2b(0,5) should throw exception" << endl;
    cout << A2b(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  vector<unsigned int> sv(3);
  sv[0] = 1;
  sv[1] = 2;
  sv[2] = 3;
  try {
    cout << "Constructor with invalid sizes should throw exception" << endl;
    RMatrix B(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  vector<int> svv(3);
  svv[0] = 1;
  svv[1] = 2;
  svv[2] = 3;
  try {
    cout << "Constructor with invalid sizes should throw exception" << endl;
    RMatrix B(svv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.clear();
  try {
    cout << "Constructor with invalid sizes should throw exception" << endl;
    RMatrix B(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.resize(1);
  sv[0] = 4;
  try {
    cout << "Constructor with invalid sizes should throw exception" << endl;
    RMatrix B(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.resize(2);
  sv[0] = 4;
  sv[1] = 0;
  {
    cout << "Constructor with a zero size should be okay" << endl;
    RMatrix B(sv);
    cout << "B size should be 0:  result = " << B.size() << endl;
    cout << "B.size(0) should be 0: result = " << B.size(0) << endl;
    cout << "B.size(1) should be 0: result = " << B.size(1) << endl;
  }

  sv.resize(2);
  sv[0] = 4;
  sv[1] = 3;
  RMatrix B0(sv);
  cout << "B0 size should be 12:  result = " << B0.size() << endl;
  cout << "B0.size(0) should be 4: result = " << B0.size(0) << endl;
  cout << "B0.size(1) should be 3: result = " << B0.size(1) << endl;

  svv.resize(2);
  svv[0] = -2;
  svv[1] = 5;
  try {
    cout << "Constructor with a negative size should throw exception" << endl;
    RMatrix B1(svv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  svv[0] = 2;
  RMatrix B1(svv);
  cout << "B1 size should be 10:  result = " << B1.size() << endl;
  cout << "B1.size(0) should be 2: result = " << B1.size(0) << endl;
  cout << "B1.size(1) should be 5: result = " << B1.size(1) << endl;

  RMatrix B2(sv, -75.3);
  for (int j = 0; j < int(B2.size(0)); ++j)
    for (int k = 0; k < int(B2.size(1)); ++k) {
      cout << "B2 should be initialized to -75.3:  result(" << j << "," << k
           << ") = " << B2(j, k) << endl;
    }

  RMatrix B3(svv, 22.3);
  for (int j = 0; j < int(B3.size(0)); ++j)
    for (int k = 0; k < int(B3.size(1)); ++k) {
      cout << "B3 should be initialized to 22.3:  result(" << j << "," << k
           << ") = " << B3(j, k) << endl;
    }

  RMatrix B4(B2);
  for (int j = 0; j < int(B4.size(0)); ++j)
    for (int k = 0; k < int(B4.size(1)); ++k) {
      cout << "B4 should be initialized to -75.3:  result(" << j << "," << k
           << ") = " << B4(j, k) << endl;
    }

  vector<unsigned int> bsizes = B4.sizes();
  cout << "B4 size(0) should be 4: result = " << bsizes[0] << endl;
  cout << "B4 size(1) should be 3: result = " << bsizes[1] << endl;

  B4.clear();
  cout << "B4 size should be zero after clear: result = " << B4.size() << endl;
  cout << "B4 size(0) should be 0: result = " << B4.size(0) << endl;
  cout << "B4 size(1) should be 0: result = " << B4.size(1) << endl;

  A2.resize();
  cout << "A2 size should be zero after clear: result = " << A2.size() << endl;
  cout << "A2 size(0) should be 0: result = " << A2.size(0) << endl;
  cout << "A2 size(1) should be 0: result = " << A2.size(1) << endl;

  A2.resize(3, 7);
  cout << "A2 size should be 21 after resize: result = " << A2.size() << endl;
  cout << "A2 size(0) should be 3: result = " << A2.size(0) << endl;
  cout << "A2 size(1) should be 7: result = " << A2.size(1) << endl;

  sv.resize(3);
  sv[0] = 1;
  sv[1] = 2;
  sv[2] = 3;
  try {
    cout << "Resize with invalid sizes should throw exception" << endl;
    A2.resize(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  svv.resize(3);
  svv[0] = 1;
  svv[1] = 2;
  svv[2] = 3;
  try {
    cout << "Resize with invalid sizes should throw exception" << endl;
    A2.resize(svv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.clear();
  try {
    cout << "Resize with invalid sizes should throw exception" << endl;
    A2.resize(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.resize(1);
  sv[0] = 4;
  try {
    cout << "Resize with invalid sizes should throw exception" << endl;
    A2.resize(sv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  sv.resize(2);
  sv[0] = 4;
  sv[1] = 0;
  {
    cout << "Constructor with a zero size should be okay" << endl;
    A2.resize(sv);
    cout << "A2 size should be 0:  result = " << A2.size() << endl;
    cout << "A2.size(0) should be 0: result = " << A2.size(0) << endl;
    cout << "A2.size(1) should be 0: result = " << A2.size(1) << endl;
  }

  sv.resize(2);
  sv[0] = 4;
  sv[1] = 3;
  A2.resize(sv);
  cout << "A2 size should be 12:  result = " << A2.size() << endl;
  cout << "A2.size(0) should be 4: result = " << A2.size(0) << endl;
  cout << "A2.size(1) should be 3: result = " << A2.size(1) << endl;

  svv.resize(2);
  svv[0] = -2;
  svv[1] = 5;
  try {
    cout << "Resize with a negative size should throw exception" << endl;
    A2.resize(svv);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  svv[0] = 2;
  A2.resize(svv);
  cout << "A2 size should be 10:  result = " << A2.size() << endl;
  cout << "A2.size(0) should be 2: result = " << A2.size(0) << endl;
  cout << "A2.size(1) should be 5: result = " << A2.size(1) << endl;

  A2(1, 4) = 46.3;
  cout << "A2(1,4) should be 46.3: result = " << A2(1, 4) << endl;
  A2.put(0, 3, 4.3);
  cout << "A2(0,3) should be 4.3: result = " << A2.get(0, 3) << endl;
  sv[0] = 1;
  sv[1] = 2;
  A2(sv) = 77.2;
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 77.2: result = " << A2(sv) << endl;
  A2.put(sv, 55.2);
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 55.2: result = " << A2.get(sv) << endl;

  ival1 = 12;
  ival2 = 9;
  A2.resize(ival1, ival2);
  svv[0] = 6;
  svv[1] = 7;
  A2(svv) = 77.2;
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 77.2: result = " << A2(svv) << endl;
  A2.put(svv, 55.2);
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 55.2: result = " << A2.get(svv) << endl;

  uival1 = 3;
  uival2 = 6;
  A2.resize(uival1, uival2);
  A2 = 25.7;
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "A2 should be set to 25.7:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 += 6.6;
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "A2 should be set to 32.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 -= 10;
  for (uint j = 0; j < A2.size(0); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "A2 should be set to 22.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 *= 10;
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (uint k = 0; k < A2.size(1); ++k) {
      cout << "A2 should be set to 223.0:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 /= 100;
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "A2 should be set to 2.230:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }

  A2.resize(3, 4);
  A1.resize(3, 4);
  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      A2(j, k) = 3.4 * (j + k);
      A1(j, k) = -2.7 * (j - 3 * k);
    }

  for (int j = 0; j < int(A2.size(0)); ++j)
    for (int k = 0; k < int(A2.size(1)); ++k) {
      cout << "A2(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A2(j, k) << endl;
      cout << "A1(" << j << "," << k << ") should be " << -2.7 * (j - 3 * k)
           << ": result = " << A1(j, k) << endl;
    }

  RMatrix A3;
  A3 = A2;
  for (int j = 0; j < int(A3.size(0)); ++j)
    for (int k = 0; k < int(A3.size(1)); ++k) {
      cout << "A3(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A3(j, k) << endl;
    }
  A3 += A1;
  for (int j = 0; j < int(A3.size(0)); ++j)
    for (int k = 0; k < int(A3.size(1)); ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) - 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }
  A3 = A2;
  A3 -= A1;
  for (int j = 0; j < int(A3.size(0)); ++j)
    for (int k = 0; k < int(A3.size(1)); ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) + 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }

  A2.resize(6, 7);
  try {
    cout << "Trying to add mismatched size: should throw exception" << endl;
    A1 += A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to subtract mismatched size: should throw exception"
         << endl;
    A1 -= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }

  uival1 = 8;
  uival2 = 6;
  A3.resize(uival1, uival2);
  int count = 0;
  for (uint k = 0; k < A3.size(1); ++k)
    for (uint j = 0; j < A3.size(0); ++j)
      A3.put(j, k, count++);

  // A3.output();  add below to Matrix class for check
  // void output()
  //{ for (uint k=0;k<m_store.size();++k) std::cout << "m_store["<<k<<"] =
  //"<<m_store[k]<<std::endl;
  // }

  cout << "simple tests done" << endl << endl;
}

template <typename T>
bool dotestMatrix(const vector<T>& datacmp, const vector<unsigned int>& idata) {
  unsigned int successcount = 0;
  unsigned int failcount = 0;
  unsigned int datacount = 0;
  unsigned int icount = 0;
  unsigned int dmax = datacmp.size();
  unsigned int imax = idata.size();

  Matrix<T>* Aptr = 0;
  vector<unsigned int> csizes;

  unsigned int ndim = 2;
  Aptr = new Matrix<T>(2, 3);
  csizes.resize(ndim);
  csizes[0] = 2;
  csizes[1] = 3;

  unsigned int sizecheck = Aptr->size();
  unsigned int ss = 1;
  for (uint k = 0; k < ndim; ++k)
    ss *= csizes[k];
  if (ss == sizecheck)
    successcount++;
  else {
    failcount++;
    cout << "size test failed" << endl;
  }
  const std::vector<uint>& scheck = Aptr->sizes();
  if (scheck == csizes)
    successcount++;
  else {
    failcount++;
    cout << "sizes test failed" << endl;
  }
  for (uint i = 0; i < 10; i++) {
    try {
      unsigned int isize = Aptr->size(i);
      if ((i < ndim) && (isize == csizes[i]))
        successcount++;
      else {
        failcount++;
        cout << "size(i) test failed" << endl;
      }
    } catch (const std::exception& xp) {
      if (i >= ndim)
        successcount++;
      else {
        failcount++;
        cout << "size(i) test failed" << endl;
      }
    }
  }
  delete Aptr;
  Aptr = 0;

  csizes.resize(ndim);
  for (uint i = 0; i < ndim; ++i)
    csizes[i] = idata[inc(icount, imax)];
  Matrix<T> A0(csizes);

  for (int jj = 0; jj < 3; jj++) {

    if (jj != 0) {
      for (uint i = 0; i < ndim; ++i)
        csizes[i] = idata[inc(icount, imax)];
      A0.resize(csizes);
    }

    unsigned int sizecheck = A0.size();
    unsigned int ss = 1;
    for (uint k = 0; k < ndim; ++k)
      ss *= csizes[k];
    if (ndim == 0)
      ss = 0;
    if (ss == sizecheck)
      successcount++;
    else {
      failcount++;
      cout << "size test failed" << endl;
    }
    const std::vector<uint>& scheck = A0.sizes();
    if (scheck == csizes)
      successcount++;
    else {
      failcount++;
      cout << "sizes test failed" << endl;
    }
    for (uint i = 0; i < 10; i++) {
      try {
        unsigned int isize = A0.size(i);
        if ((i < ndim) && (isize == csizes[i]))
          successcount++;
        else {
          failcount++;
          cout << "size(i) test failed" << endl;
        }
      } catch (const std::exception& xp) {
        if (i >= ndim)
          successcount++;
        else {
          failcount++;
          cout << "size(i) test failed" << endl;
        }
      }
    }

    datacount = 0;
    unsigned int dmax = datacmp.size();
    for (uint k0 = 0; k0 < csizes[0]; ++k0)
      for (uint k1 = 0; k1 < csizes[1]; ++k1) {
        A0(k0, k1) = datacmp[inc(datacount, dmax)];
        if (datacount >= dmax)
          datacount = 0;
      }

    datacount = 0;
    for (uint k0 = 0; k0 < csizes[0]; ++k0)
      for (uint k1 = 0; k1 < csizes[1]; ++k1) {
        if (A0(k0, k1) != datacmp[inc(datacount, dmax)])
          failcount++;
        else
          successcount++;
        if (datacount >= dmax)
          datacount = 0;
      }

    datacount = 0;
    dmax = datacmp.size();
    for (uint k0 = 0; k0 < csizes[0]; ++k0)
      for (uint k1 = 0; k1 < csizes[1]; ++k1) {
        A0.put(k0, k1, datacmp[inc(datacount, dmax)]);
        if (datacount >= dmax)
          datacount = 0;
      }

    datacount = 0;
    for (uint k0 = 0; k0 < csizes[0]; ++k0)
      for (uint k1 = 0; k1 < csizes[1]; ++k1) {
        if (A0.get(k0, k1) != datacmp[inc(datacount, dmax)])
          failcount++;
        else
          successcount++;
        if (datacount >= dmax)
          datacount = 0;
      }
  }

  vector<unsigned int> isizes(ndim);
  for (uint k = 0; k < ndim; ++k)
    isizes[k] = idata[inc(icount, imax)];
  T cf = datacmp[inc(datacount, dmax)];
  Matrix<T> B0(isizes, cf);
  MultiIntLooper ind(isizes);
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  cf = datacmp[inc(datacount, dmax)];
  cf *= datacmp[inc(datacount, dmax)];
  cf -= datacmp[inc(datacount, dmax)];
  B0 = cf;
  cout << "cf = " << cf << endl;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  T cf2 = datacmp[inc(datacount, dmax)];
  cf += cf2;
  B0 += cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf *= cf2;
  B0 *= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf -= cf2;
  B0 -= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  cf2 = datacmp[inc(datacount, dmax)];
  cf /= cf2;
  B0 /= cf2;
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(B0(ind.get_current()), cf))
      successcount++;
    else
      failcount++;
  }

  Matrix<T> W(B0);
  for (ind.start(); ind.notdone(); ++ind) {
    if (equal(W(ind.get_current()), B0(ind.get_current())))
      successcount++;
    else
      failcount++;
  }

  for (uint k = 0; k < ndim; ++k)
    isizes[k] = idata[inc(icount, imax)];
  B0.resize(isizes);
  MultiIntLooper dil(isizes);
  datacount = 0;
  for (dil.start(); dil.notdone(); ++dil)
    B0(dil.get_current()) = datacmp[inc(datacount, dmax)];
  datacount = dmax / 2;
  Matrix<T> B1(isizes);
  for (dil.start(); dil.notdone(); ++dil)
    B1(dil.get_current()) = datacmp[inc(datacount, dmax)];

  Matrix<T> B2(B0);
  B2 += B1;
  datacount = 0;
  unsigned int dcount2 = dmax / 2;
  for (dil.start(); dil.notdone(); ++dil) {
    if (equal(B2(dil.get_current()),
              (datacmp[inc(datacount, dmax)] + datacmp[inc(dcount2, dmax)])))
      successcount++;
    else
      failcount++;
  }

  cout << endl;
  cout << " Number of successful tests = " << successcount << endl;
  cout << "     Number of FAILED tests = " << failcount << endl;
  return (failcount == 0);
}

void testMatrix() {

  cout << "Starting Matrix tests" << endl;

  dotestMatrixSimple();

  bool success = true;
  srand(time(NULL));
  unsigned int nintegers = 2056;
  vector<unsigned int> idata(nintegers);
  for (unsigned int k = 0; k < nintegers; ++k)
    idata[k] = rand() % 8 + 1;

  unsigned int ndata = 5048;
  vector<int> intdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    intdata[k] = (rand() % 4096 + 1) - 2047;

  int ndim = 2;
  cout << endl << "Tests: integer type, ndim = " << ndim << endl;
  success = success && dotestMatrix(intdata, idata);
  cout << endl;

  ndata = 5048;
  vector<unsigned int> uintdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    uintdata[k] = (rand() % 4096 + 1);

  cout << endl << "Tests: unsigned integer type, ndim = " << ndim << endl;
  success = success && dotestMatrix(uintdata, idata);
  cout << endl;

  ndata = 5048;
  vector<float> fdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    fdata[k] =
        1.434 * ((rand() % 4096 + 1) - 2047) - 0.54511 * ((rand() % 128) - 63);

  cout << endl << "Tests: float type, ndim = " << ndim << endl;
  success = success && dotestMatrix(fdata, idata);
  cout << endl;

  ndata = 5048;
  vector<double> ddata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    ddata[k] = -2.3418 * ((rand() % 4096 + 1) - 2047) -
               0.22485 * ((rand() % 128) - 63);

  cout << endl << "Tests: double type, ndim = " << ndim << endl;
  success = success && dotestMatrix(ddata, idata);
  cout << endl;

  ndata = 5048;
  vector<complex<float>> zfdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    zfdata[k] = complex<float>(1.434 * ((rand() % 4096 + 1) - 2047),
                               -0.54511 * ((rand() % 128) - 63));

  cout << endl << "Tests: complex<float> type, ndim = " << ndim << endl;
  success = success && dotestMatrix(zfdata, idata);
  cout << endl;

  ndata = 5048;
  vector<complex<double>> zddata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    zddata[k] = complex<double>(-2.3418 * ((rand() % 4096 + 1) - 2047),
                                0.22485 * ((rand() % 128) - 63));

  cout << endl << "Tests: complex<double> type, ndim = " << ndim << endl;
  success = success && dotestMatrix(zddata, idata);
  cout << endl;

  CVector vt(5), vs(5);
  vt[0] = complex<double>(1.43, -5.234);
  vt[1] = complex<double>(8.321, 0.432);
  vt[2] = complex<double>(12.324, -42.11);
  vt[3] = complex<double>(-6.32, 24.11);
  vt[4] = complex<double>(0.543, 12.234);
  vs[0] = complex<double>(1.324, 9.546);
  vs[1] = complex<double>(5.23, 8.32);
  vs[2] = complex<double>(-4.442, 3.33);
  vs[3] = complex<double>(0.412, -3.427);
  vs[4] = complex<double>(7.34, 22.22);

  complex<double> z = dotProduct(vs, vt);
  z = z - complex<double>(-5.330592000, 124.4685240);
  if (std::abs(z) > 1e-12)
    success = false;

  ComplexHermitianMatrix At(5);
  At.put(0, 0, complex<double>(18.0, 0.0));
  At.put(0, 1, complex<double>(56.0, -79.0));
  At.put(0, 2, complex<double>(-40.0, 161.0));
  At.put(0, 3, complex<double>(-35.0, 84.0));
  At.put(0, 4, complex<double>(76.0, 154.0));
  At.put(1, 1, complex<double>(188.0, 0.0));
  At.put(1, 2, complex<double>(-4.0, 46.0));
  At.put(1, 3, complex<double>(-83.0, 60.0));
  At.put(1, 4, complex<double>(-48.0, 11.0));
  At.put(2, 2, complex<double>(-18.0, 0.0));
  At.put(2, 3, complex<double>(-88.0, -24.0));
  At.put(2, 4, complex<double>(-99.0, -36.0));
  At.put(3, 3, complex<double>(-36.0, 0.0));
  At.put(3, 4, complex<double>(144.0, -113.0));
  At.put(4, 4, complex<double>(54.0, 0.0));

  cout.precision(14);
  CVector outvec;
  multiply(outvec, At, vt);
  z = outvec[0] - complex<double>(3165.7860, 2579.8610);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[1] - complex<double>(2863.0, -2325.1630);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[2] - complex<double>(386.3490, -2848.09800);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[3] - complex<double>(459.8530, 4361.78700);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[4] - complex<double>(-4401.316, 7300.59900);
  if (std::abs(z) > 1e-12)
    success = false;

  z = dotProduct(vt, outvec);
  z = z - complex<double>(327719.433134, 0.0);
  if (std::abs(z) > 1e-6)
    success = false;
  z = dotProduct(vs, outvec);
  z = z - complex<double>(128401.924739, 103335.700119);
  if (std::abs(z) > 1e-6)
    success = false;

  multiply(outvec, At, vs);
  z = outvec[0] - complex<double>(-1975.0500, 2349.84900);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[1] - complex<double>(-257.47800, 1309.02100);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[2] - complex<double>(1880.45800, -3101.136000);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[3] - complex<double>(4684.59800, 644.2980000);
  if (std::abs(z) > 1e-12)
    success = false;
  z = outvec[4] - complex<double>(2574.00500, 328.0760000);
  if (std::abs(z) > 1e-12)
    success = false;

  RealSymmetricMatrix Ar(5);
  Ar(0, 0) = 3.28358208947999985;
  Ar(0, 1) = -1.71641791040999991;
  Ar(0, 2) = 9.92537313411000000;
  Ar(0, 3) = 1.11940298505000002;
  Ar(0, 4) = 3.50746268648999981;
  Ar(1, 1) = 3.13432835813999988;
  Ar(1, 2) = 0.597014925360000004;
  Ar(1, 3) = 11.1940298504999998;
  Ar(1, 4) = 7.53731343267000042;
  Ar(2, 2) = -10.4477611937999999;
  Ar(2, 3) = 0.671641791029999990;
  Ar(2, 4) = 12.1641791042099996;
  Ar(3, 3) = -0.14925373134000000;
  Ar(3, 4) = 7.01492537297999963;
  Ar(4, 4) = 2.68656716411999996;

  RVector rt(5);
  rt[0] = -2.65895953757225433;
  rt[1] = -7.28323699421965317;
  rt[2] = -3.00578034682080925;
  rt[3] = 3.46820809248554913;
  rt[4] = 1.15606936416184971;

  RVector rs(5);
  rs[0] = -7.26577437858508603;
  rs[1] = -7.26577437858508603;
  rs[2] = 17.3996175908221797;
  rs[3] = -0.19120458891013384;
  rs[4] = 12.0458891013384321;

  RVector ovec;
  multiply(ovec, Ar, rt);
  double rz;
  rz = ovec[0] + 18.1261323436612721;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[1] - 27.4782158565260119;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[2] - 17.0563368126693634;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[3] + 78.9319299439109816;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[4] + 73.3500129394612741;
  if (std::abs(rz) > 1e-12)
    success = false;

  rz = dotProduct(rt, ovec);
  rz = rz + 561.751368774367065;
  if (std::abs(rz) > 1e-9)
    success = false;
  rz = dotProduct(rs, ovec);
  rz = rz + 639.650364967131420;
  if (std::abs(rz) > 1e-9)
    success = false;

  multiply(ovec, Ar, rs);
  rz = ovec[0] - 203.347507201370935;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[1] - 88.7389058511912098;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[2] + 111.840415510795414;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[3] - 6.7492366083957903;
  if (std::abs(rz) > 1e-12)
    success = false;
  rz = ovec[4] - 162.424017575833642;
  if (std::abs(rz) > 1e-12)
    success = false;

  if (success)
    cout << "ALL MATRIX TESTS PASSED!!" << endl;
  else
    cout << "Some Matrix tests FAILED" << endl;
}

// ****************************************************************************

void dotestComplexHermitianMatrixSimple() {
  cout << endl << endl << "Simple Test Results" << endl << endl;

  ComplexHermitianMatrix A0;
  cout << "A0 size should be 0:  result = " << A0.size() << endl;
  try {
    cout << "Trying A0(2,3): should throw exception" << endl;
    cout << A0(2, 3) << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  int ival1 = 4, ival2 = 3;
  ComplexHermitianMatrix A1(ival1);
  cout << "A1 size should be 4:  result = " << A1.size() << endl;
  uint uival1 = 17;
  ComplexHermitianMatrix A1b(uival1);
  cout << "A1b size should be 17:  result = " << A1b.size() << endl;

  ival1 = 3;
  ival2 = 2;
  ComplexHermitianMatrix A2(ival1, 4.6);
  cout << "A2 size should be 3:  result = " << A2.size() << endl;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "Should be initialized to 4.6:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  try {
    cout << "A2(4,5) should throw exception" << endl;
    cout << A2(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2(0,5) should throw exception" << endl;
    cout << A2(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  uival1 = 5;
  ComplexHermitianMatrix A2b(uival1, 24.6);
  cout << "A2b size should be 5:  result = " << A2b.size() << endl;
  for (uint j = 0; j < A2b.size(); ++j)
    for (uint k = 0; k < A2b.size(); ++k) {
      cout << "Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  try {
    cout << "A2b(4,5) should throw exception" << endl;
    cout << A2b(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2b(0,5) should throw exception" << endl;
    cout << A2b(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  int count = 0;
  A2.resize(8);
  for (int k = 0; k < int(A2.size()); ++k)
    for (int j = 0; j <= k; ++j) {
      A2.put(j, k, complex<double>(double(count), double(6.7123 * (k - j))));
      count++;
    }
  ival1 = 2;
  ival2 = 7;
  A2.put(ival2, ival1, complex<double>(456.25, 33.2));
  cout << "A2(2,7) should be (456.25,-33.2): result = " << A2(ival1, ival2)
       << endl;
  cout << "A2(7,2) should be (456.25,33.2): result = " << A2.get(ival2, ival1)
       << endl;
  uival1 = 3;
  int uival2 = 5;
  A2.put(uival1, uival2, complex<double>(723.9, -89.32));
  cout << "A2(3,5) should be (723.9,-89.32): result = "
       << A2.get(uival1, uival2) << endl;
  cout << "A2(5,3) should be (723.9,89.32): result = " << A2(uival2, uival1)
       << endl;
  uival1 = 2;
  try {
    cout << "attempt to assign nonzero imaginary part to diagonal element: "
            "should throw exception"
         << endl;
    A2.put(uival1, uival1, complex<double>(723.9, -89.32));
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "exception caught: correct" << endl;
  }
  A2.put(uival1, uival1, complex<double>(62.7, 0.0));
  ival1 = 2;
  cout << "A2(2,2) should be (62.7, 0): result = " << A2.get(uival1, uival1)
       << endl;
  cout << "A2(2,2) should be (62.7, 0): result = " << A2(ival1, ival1) << endl;

  // A2.output();  //add below to Matrix class for check
  // void output()
  //{ for (uint k=0;k<m_store.size();++k) std::cout << "m_store["<<k<<"] =
  //"<<m_store[k]<<std::endl;
  //}

  try {
    cout << "Constructor with a negative size should throw exception" << endl;
    ComplexHermitianMatrix B1(-3);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  ComplexHermitianMatrix B1(0);
  cout << "B1 size should be 0:  result = " << B1.size() << endl;

  ComplexHermitianMatrix A4(A2);
  for (int j = 0; j < int(A4.size()); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }

  A4.clear();
  cout << "A4 size should be zero after clear: result = " << A4.size() << endl;

  A4 = A2;
  for (int j = 0; j < int(A4.size()); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (uint j = 0; j < A4.size(); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (int j = 0; j < int(A4.size()); ++j)
    for (uint k = 0; k < A4.size(); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (uint j = 0; j < A4.size(); ++j)
    for (uint k = 0; k < A4.size(); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }

  try {
    cout << "resize with negative value should throw exception" << endl;
    A2.resize(-7);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  A2.resize(5);
  cout << "A2 resize should be 5:  result = " << A2.size() << endl;

  A2.put(1, 4, 46.3);
  cout << "A2(1,4) should be 46.3: result = " << A2(1, 4) << endl;
  A2.put(0, 3, 4.3);
  cout << "A2(0,3) should be 4.3: result = " << A2.get(0, 3) << endl;
  vector<int> sv(3);
  sv[0] = 1;
  sv[1] = 2;
  sv[2] = 2;
  try {
    cout << "invalid element access" << endl;
    cout << A2(sv) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  sv.resize(2);
  sv[0] = 1;
  sv[1] = 2;
  A2.put(sv, 77.2);
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 77.2: result = " << A2(sv) << endl;
  A2.put(sv, 55.2);
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 55.2: result = " << A2.get(sv) << endl;

  ival2 = 9;
  A2.resize(ival2);
  vector<unsigned int> svv(1);
  svv[0] = 5;
  try {
    cout << "invalid element access" << endl;
    cout << A2(svv) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  svv.resize(2);
  svv[0] = 6;
  svv[1] = 7;
  A2.put(svv, 77.2);
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 77.2: result = " << A2(svv) << endl;
  A2.put(svv, 55.2);
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 55.2: result = " << A2.get(svv) << endl;

  uival2 = 6;
  A2.resize(uival2);
  A2 = 25.7;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 25.7:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 += 6.6;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 32.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 -= 10;
  for (uint j = 0; j < A2.size(); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 22.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 *= 10;
  for (int j = 0; j < int(A2.size()); ++j)
    for (uint k = 0; k < A2.size(); ++k) {
      cout << "A2 should be set to 223.0:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 /= 100;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 2.230:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }

  A2.resize(4);
  A1.resize(4);
  for (int k = 0; k < int(A2.size()); ++k)
    for (int j = 0; j <= k; ++j) {
      A2.put(j, k, complex<double>(3.4 * (j + k), -1.2 * (j - k)));
      A1.put(j, k, complex<double>(-2.7 * (j + k), -5.2 * (k - j)));
    }

  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2(" << j << "," << k << ") should be (" << 3.4 * (j + k) << ","
           << -1.2 * (j - k) << "): result = " << A2(j, k) << endl;
      cout << "A1(" << j << "," << k << ") should be (" << -2.7 * (j + k) << ","
           << -5.2 * (k - j) << "): result = " << A1(j, k) << endl;
    }

  ComplexHermitianMatrix A3;
  A3 = A2;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "(1) A3(" << j << "," << k << ") should be (" << 3.4 * (j + k)
           << "," << -1.2 * (j - k) << "): result = " << A3(j, k) << endl;
    }
  A3 += A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "(2) A3(" << j << "," << k << ") should be ("
           << 3.4 * (j + k) - 2.7 * (j + k) << ","
           << -1.2 * (j - k) - 5.2 * (k - j) << "): result = " << A3(j, k)
           << endl;
    }
  A3 = A2;
  A3 -= A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "(3) A3(" << j << "," << k << ") should be ("
           << 3.4 * (j + k) + 2.7 * (j + k) << ","
           << -1.2 * (j - k) + 5.2 * (k - j) << "): result = " << A3(j, k)
           << endl;
    }

  A2.resize(7);
  try {
    cout << "Trying to add mismatched size: should throw exception" << endl;
    A1 += A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to subtract mismatched size: should throw exception"
         << endl;
    A1 -= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }

  cout << "simple tests done" << endl << endl;
}

bool dotestComplexHermitianMatrix(const vector<complex<double>>& datacmp,
                                  const vector<unsigned int>& idata) {
  unsigned int successcount = 0;
  unsigned int failcount = 0;

  int firstsize = 6;
  ComplexHermitianMatrix H0(firstsize);
  for (int row = 0; row < firstsize; row++)
    for (int col = row; col < firstsize; col++)
      H0.put(row, col, 100 * row + col);

  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(double(100 * row + col), H0(row, col)))
        successcount++;
      else
        failcount++;

  // cout << "cf="<<datacmp[0]<<endl;
  ComplexHermitianMatrix H1(firstsize, real(datacmp[0]));
  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(H1(row, col), real(datacmp[0])))
        successcount++;
      else
        failcount++;

  unsigned int secondsize =
      idata[3] * 128; // cout << "secondsize = "<<secondsize<<endl;
  ComplexHermitianMatrix H2(secondsize);
  unsigned int count = 0;
  for (int row = 0; row < int(secondsize); row++)
    for (int col = row; col < int(secondsize); col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      H2.put(row, col, temp);
      if (count > datacmp.size())
        count = 0;
    }

  count = 0;
  for (int row = 0; row < int(secondsize); row++)
    for (int col = row; col < int(secondsize); col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H2(row, col), temp))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  ComplexHermitianMatrix H3;
  H3 = H2;
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H3(row, col), temp))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  ComplexHermitianMatrix H4(H2);
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H3(row, col), temp))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H1;
  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(H4(row, col), real(datacmp[0])))
        successcount++;
      else
        failcount++;

  H4.resize(4);
  H4 = imag(datacmp[5]);
  for (int row = 0; row < 4; row++)
    for (int col = row; col < 4; col++) {
      if (equal(H4(row, col), imag(datacmp[5])))
        successcount++;
      else
        failcount++;
    }

  H4.clear();
  if (H4.size() == 0)
    successcount++;
  else
    failcount++;

  H3 += real(datacmp[8]);
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H3(row, col), temp + real(datacmp[8])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 -= real(datacmp[6]);
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H3(row, col), temp + real(datacmp[8] - datacmp[6])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 *= imag(datacmp[2]);
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H3(row, col),
                (temp + real(datacmp[8] - datacmp[6])) * imag(datacmp[2])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H3;
  H4 += H3;
  count = 0;
  for (uint row = 0; row < secondsize; row++)
    for (uint col = row; col < secondsize; col++) {
      complex<double> temp(datacmp[count++]);
      if (row == col)
        temp = complex<double>(real(temp), 0.0);
      if (equal(H4(row, col), (temp + real(datacmp[8] - datacmp[6])) * 2.0 *
                                  imag(datacmp[2])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  cout << endl;
  cout << " Number of successful tests = " << successcount << endl;
  cout << "     Number of FAILED tests = " << failcount << endl;
  return (failcount == 0);
}

void testComplexHermitianMatrix() {
  bool success = true;

  dotestComplexHermitianMatrixSimple();

  srand(time(NULL));
  unsigned int nintegers = 2056;
  vector<unsigned int> idata(nintegers);
  for (unsigned int k = 0; k < nintegers; ++k)
    idata[k] = rand() % 8 + 1;

  unsigned int ndata = 5048;
  vector<complex<double>> fdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    fdata[k] = complex<double>(1.434 * ((rand() % 4096 + 1) - 2047) -
                                   0.54511 * ((rand() % 128) - 63),
                               3.781 * ((rand() % 4096 + 1) - 2047));

  success = dotestComplexHermitianMatrix(fdata, idata);

  if (success)
    cout << "ALL ComplexHermitianMatrix TESTS PASSED!!" << endl;
  else
    cout << "Some ComplexHermitianMatrix tests FAILED" << endl;
}

void dotestRealSymmetricMatrixSimple() {
  cout << endl << endl << "Simple Test Results" << endl << endl;

  RealSymmetricMatrix A0;
  cout << "A0 size should be 0:  result = " << A0.size() << endl;
  try {
    cout << "Trying A0(2,3): should throw exception" << endl;
    cout << A0(2, 3) << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  int ival1 = 4, ival2 = 3;
  RealSymmetricMatrix A1(ival1);
  cout << "A1 size should be 4:  result = " << A1.size() << endl;
  uint uival1 = 17;
  RealSymmetricMatrix A1b(uival1);
  cout << "A1b size should be 17:  result = " << A1b.size() << endl;

  ival1 = 3;
  ival2 = 2;
  RealSymmetricMatrix A2(ival1, 4.6);
  cout << "A2 size should be 3:  result = " << A2.size() << endl;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "Should be initialized to 4.6:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  try {
    cout << "A2(4,5) should throw exception" << endl;
    cout << A2(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2(0,5) should throw exception" << endl;
    cout << A2(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  uival1 = 5;
  RealSymmetricMatrix A2b(uival1, 24.6);
  cout << "A2b size should be 5:  result = " << A2b.size() << endl;
  for (uint j = 0; j < A2b.size(); ++j)
    for (uint k = 0; k < A2b.size(); ++k) {
      cout << "Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  try {
    cout << "A2b(4,5) should throw exception" << endl;
    cout << A2b(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2b(0,5) should throw exception" << endl;
    cout << A2b(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  int count = 0;
  A2.resize(8);
  for (int k = 0; k < int(A2.size()); ++k)
    for (int j = 0; j <= k; ++j) {
      A2(j, k) = count++;
    }
  ival1 = 2;
  ival2 = 7;
  A2(ival2, ival1) = 456.25;
  cout << "A2(2,7) should be 456.25: result = " << A2(ival1, ival2) << endl;
  cout << "A2(2,7) should be 456.25: result = " << A2.get(ival1, ival2) << endl;
  uival1 = 3;
  int uival2 = 5;
  A2.put(uival2, uival1, 723.9);
  cout << "A2(5,3) should be 723.9: result = " << A2.get(uival1, uival2)
       << endl;
  cout << "A2(5,3) should be 723.9: result = " << A2(uival1, uival2) << endl;

  // A2.output();  //add below to Matrix class for check
  // void output()
  //{ for (uint k=0;k<m_store.size();++k) std::cout << "m_store["<<k<<"] =
  //"<<m_store[k]<<std::endl;
  //}

  try {
    cout << "Constructor with a negative size should throw exception" << endl;
    RealSymmetricMatrix B1(-3);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  RealSymmetricMatrix B1(0);
  cout << "B1 size should be 0:  result = " << B1.size() << endl;

  RealSymmetricMatrix A4(A2);
  for (int j = 0; j < int(A4.size()); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }

  A4.clear();
  cout << "A4 size should be zero after clear: result = " << A4.size() << endl;

  A4 = A2;
  for (int j = 0; j < int(A4.size()); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (uint j = 0; j < A4.size(); ++j)
    for (int k = 0; k < int(A4.size()); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (int j = 0; j < int(A4.size()); ++j)
    for (uint k = 0; k < A4.size(); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }
  for (uint j = 0; j < A4.size(); ++j)
    for (uint k = 0; k < A4.size(); ++k) {
      cout << "A4(" << j << "," << k << ") should be initialized be "
           << A2(j, k) << ":  result = " << A4(j, k) << endl;
    }

  try {
    cout << "resize with negative value should throw exception" << endl;
    A2.resize(-7);
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  A2.resize(5);
  cout << "A2 resize should be 5:  result = " << A2.size() << endl;

  A2(1, 4) = 46.3;
  cout << "A2(1,4) should be 46.3: result = " << A2(1, 4) << endl;
  A2.put(0, 3, 4.3);
  cout << "A2(0,3) should be 4.3: result = " << A2.get(0, 3) << endl;
  vector<int> sv(3);
  sv[0] = 1;
  sv[1] = 2;
  sv[2] = 2;
  try {
    cout << "invalid element access" << endl;
    cout << A2(sv) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  sv.resize(2);
  sv[0] = 1;
  sv[1] = 2;
  A2(sv) = 77.2;
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 77.2: result = " << A2(sv) << endl;
  A2.put(sv, 55.2);
  cout << "A2[" << sv[0] << "," << sv[1]
       << "] should be 55.2: result = " << A2.get(sv) << endl;

  ival1 = 12;
  ival2 = 9;
  A2.resize(ival2);
  vector<unsigned int> svv(1);
  svv[0] = 5;
  try {
    cout << "invalid element access" << endl;
    cout << A2(svv) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "caught exception: correct" << endl;
  }

  svv.resize(2);
  svv[0] = 6;
  svv[1] = 7;
  A2(svv) = 77.2;
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 77.2: result = " << A2(svv) << endl;
  A2.put(svv, 55.2);
  cout << "A2[" << svv[0] << "," << svv[1]
       << "] should be 55.2: result = " << A2.get(svv) << endl;

  uival2 = 6;
  A2.resize(uival2);
  A2 = 25.7;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 25.7:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 += 6.6;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 32.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 -= 10;
  for (uint j = 0; j < A2.size(); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 22.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 *= 10;
  for (int j = 0; j < int(A2.size()); ++j)
    for (uint k = 0; k < A2.size(); ++k) {
      cout << "A2 should be set to 223.0:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 /= 100;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2 should be set to 2.230:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }

  A2.resize(4);
  A1.resize(4);
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      A2(j, k) = 3.4 * (j + k);
      A1(j, k) = -2.7 * (j - 3 * k);
    }

  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k < int(A2.size()); ++k) {
      cout << "A2(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A2(j, k) << endl;
      cout << "A1(" << j << "," << k << ") should be " << -2.7 * (j - 3 * k)
           << ": result = " << A1(j, k) << endl;
    }

  RealSymmetricMatrix A3;
  A3 = A2;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "A3(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A3(j, k) << endl;
    }
  A3 += A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) - 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }
  A3 = A2;
  A3 -= A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k < int(A3.size()); ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) + 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }

  A2.resize(7);
  try {
    cout << "Trying to add mismatched size: should throw exception" << endl;
    A1 += A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to subtract mismatched size: should throw exception"
         << endl;
    A1 -= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }

  cout << "simple tests done" << endl << endl;
}

bool dotestRealSymmetricMatrix(const vector<double>& datacmp,
                               const vector<unsigned int>& idata) {
  unsigned int successcount = 0;
  unsigned int failcount = 0;

  int firstsize = 6;
  RealSymmetricMatrix H0(firstsize);
  for (int row = 0; row < firstsize; row++)
    for (int col = row; col < firstsize; col++)
      H0(row, col) = 100 * row + col;

  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(double(100 * row + col), H0(row, col)))
        successcount++;
      else
        failcount++;

  // cout << "cf="<<datacmp[0]<<endl;
  RealSymmetricMatrix H1(firstsize, datacmp[0]);
  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(H1(row, col), datacmp[0]))
        successcount++;
      else
        failcount++;

  int secondsize = idata[3] * 128;
  cout << "secondsize = " << secondsize << endl;
  RealSymmetricMatrix H2(secondsize);
  unsigned int count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      H2(row, col) = datacmp[count++];
      if (count > datacmp.size())
        count = 0;
    }

  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H2(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  RealSymmetricMatrix H3;
  H3 = H2;
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H3(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  RealSymmetricMatrix H4(H2);
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H3(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H1;
  for (int col = firstsize - 1; col >= 0; col--)
    for (int row = 0; row <= col; row++)
      if (equal(H4(row, col), datacmp[0]))
        successcount++;
      else
        failcount++;

  H4.resize(4);
  H4 = datacmp[5];
  for (int row = 0; row < 4; row++)
    for (int col = row; col < 4; col++) {
      if (equal(H4(row, col), datacmp[5]))
        successcount++;
      else
        failcount++;
    }

  H4.clear();
  if (H4.size() == 0)
    successcount++;
  else
    failcount++;

  H3 += datacmp[8];
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H3(row, col), datacmp[count++] + datacmp[8]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 -= datacmp[6];
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H3(row, col), datacmp[count++] + datacmp[8] - datacmp[6]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 *= datacmp[2];
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H3(row, col),
                datacmp[2] * (datacmp[count++] + datacmp[8] - datacmp[6])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H3;
  H4 += H3;
  count = 0;
  for (int row = 0; row < secondsize; row++)
    for (int col = row; col < secondsize; col++) {
      if (equal(H4(row, col),
                2 * datacmp[2] * (datacmp[count++] + datacmp[8] - datacmp[6])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  cout << endl;
  cout << " Number of successful tests = " << successcount << endl;
  cout << "     Number of FAILED tests = " << failcount << endl;
  return (failcount == 0);
}

void testRealSymmetricMatrix() {
  bool success = true;

  dotestRealSymmetricMatrixSimple();

  srand(time(NULL));
  unsigned int nintegers = 2056;
  vector<unsigned int> idata(nintegers);
  for (unsigned int k = 0; k < nintegers; ++k)
    idata[k] = rand() % 8 + 1;

  unsigned int ndata = 5048;
  vector<double> fdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    fdata[k] =
        1.434 * ((rand() % 4096 + 1) - 2047) - 0.54511 * ((rand() % 128) - 63);

  success = dotestRealSymmetricMatrix(fdata, idata);

  if (success)
    cout << "ALL RealSymmetricMatrix TESTS PASSED!!" << endl;
  else
    cout << "Some RealSymmetricMatrix tests FAILED" << endl;
}

// *********************************************************************************

void dotestLTMatrixSimple() {
  cout << endl << endl << "Simple Test Results" << endl << endl;

  LTMatrix A0;
  cout << "A0 size should be 0:  result = " << A0.size() << endl;
  try {
    cout << "Trying A0(2,3): should throw exception" << endl;
    cout << A0(2, 3) << endl;
    cout << "NOT caught: ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  int ival1 = 4;
  LTMatrix A1(ival1);
  cout << "A1 size should be 4:  result = " << A1.size() << endl;

  uint uival1 = 6;
  LTMatrix A1b(uival1);
  cout << "A1b size should be 6:  result = " << A1b.size() << endl;

  ival1 = 5;
  LTMatrix A2(ival1, 4.6);
  cout << "A2 size should be 5:  result = " << A2.size() << endl;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "Should be initialized to 4.6:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  try {
    cout << "A2(4,5) should throw exception" << endl;
    cout << A2(4, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2(0,5) should throw exception" << endl;
    cout << A2(0, 5) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }
  try {
    cout << "A2(2,4) should throw exception" << endl;
    cout << A2(2, 4) << endl;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught:  correct" << endl;
  }

  uival1 = 15;
  LTMatrix A2b(uival1, 24.6);
  cout << "A2b size should be 15:  result = " << A2b.size() << endl;
  for (int j = 0; j < int(A2b.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2b Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  for (uint j = 0; j < A2b.size(); ++j)
    for (int k = 0; k <= int(j); ++k) {
      cout << "A2b Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  for (uint j = 0; j < A2b.size(); ++j)
    for (uint k = 0; k <= j; ++k) {
      cout << "A2b Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }
  for (int j = 0; j < int(A2b.size()); ++j)
    for (uint k = 0; k <= uint(j); ++k) {
      cout << "A2b Should be initialized to 24.6:  result(" << j << "," << k
           << ") = " << A2b(j, k) << endl;
    }

  int count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      A2(j, k) = count;
      cout << "Should be initialized to " << count << ":  result(" << j << ","
           << k << ") = " << A2(j, k) << endl;
      count++;
    }

  // A2.output();
  // void output()
  //{ for (uint k=0;k<m_store.size();++k) std::cout << "m_store["<<k<<"] =
  //"<<m_store[k]<<std::endl;
  // }

  count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "Should be initialized to " << count << ":  result(" << j << ","
           << k << ") = " << A2.get(j, k) << endl;
      count++;
    }
  count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      A2.put(j, k, 10 * count);
      cout << "Should be initialized to " << 10 * count << ":  result(" << j
           << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }

  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (uint k = 0; k <= j; ++k) {
      A2(j, k) = 2.2 * double(count);
      cout << "Should be initialized to " << 2.2 * double(count) << ":  result("
           << j << "," << k << ") = " << A2(j, k) << endl;
      count++;
    }
  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (uint k = 0; k <= j; ++k) {
      cout << "Should be initialized to " << 2.2 * double(count) << ":  result("
           << j << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }
  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (uint k = 0; k <= j; ++k) {
      A2.put(j, k, 8.3 * double(count));
      cout << "Should be initialized to " << 8.3 * double(count) << ":  result("
           << j << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }

  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (int k = 0; k <= int(j); ++k) {
      A2(j, k) = 2.4 * double(count);
      cout << "Should be initialized to " << 2.4 * double(count) << ":  result("
           << j << "," << k << ") = " << A2(j, k) << endl;
      count++;
    }
  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (int k = 0; k <= int(j); ++k) {
      cout << "Should be initialized to " << 2.4 * double(count) << ":  result("
           << j << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }
  count = 0;
  for (uint j = 0; j < A2.size(); ++j)
    for (int k = 0; k <= int(j); ++k) {
      A2.put(j, k, 18.3 * double(count));
      cout << "Should be initialized to " << 18.3 * double(count)
           << ":  result(" << j << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }

  count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (uint k = 0; k <= uint(j); ++k) {
      A2(j, k) = count;
      cout << "Should be initialized to " << count << ":  result(" << j << ","
           << k << ") = " << A2(j, k) << endl;
      count++;
    }
  count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (uint k = 0; k <= uint(j); ++k) {
      cout << "Should be initialized to " << count << ":  result(" << j << ","
           << k << ") = " << A2.get(j, k) << endl;
      count++;
    }
  count = 0;
  for (int j = 0; j < int(A2.size()); ++j)
    for (uint k = 0; k <= uint(j); ++k) {
      A2.put(j, k, 10 * count);
      cout << "Should be initialized to " << 10 * count << ":  result(" << j
           << "," << k << ") = " << A2.get(j, k) << endl;
      count++;
    }

  A2.clear();
  cout << "A2 size should be zero after clear: result = " << A2.size() << endl;
  A2.resize(7);
  cout << "A2 size should be 7 after resize: result = " << A2.size() << endl;
  A2.resize();
  cout << "A2 size should be 0 after resize: result = " << A2.size() << endl;

  count = 0;
  uival1 = 5;
  A2.resize(uival1);
  vector<int> ind(2);
  for (ind[0] = 0; ind[0] < int(A2.size()); ++ind[0])
    for (ind[1] = 0; ind[1] <= ind[0]; ++ind[1]) {
      A2.put(ind, 12 * count);
      cout << "Should be initialized to " << 12 * count << ":  result("
           << ind[0] << "," << ind[1] << ") = " << A2.get(ind) << endl;
      count++;
    }
  count = 0;
  for (ind[0] = 0; ind[0] < int(A2.size()); ++ind[0])
    for (ind[1] = 0; ind[1] <= ind[0]; ++ind[1]) {
      A2(ind) = 7 * count;
      cout << "Should be initialized to " << 7 * count << ":  result(" << ind[0]
           << "," << ind[1] << ") = " << A2(ind) << endl;
      count++;
    }
  count = 0;
  vector<unsigned int> iind(2);
  for (iind[0] = 0; iind[0] < A2.size(); ++iind[0])
    for (iind[1] = 0; iind[1] <= iind[0]; ++iind[1]) {
      A2(iind) = 8 * count;
      cout << "Should be initialized to " << 8 * count << ":  result("
           << iind[0] << "," << iind[1] << ") = " << A2(iind) << endl;
      count++;
    }
  for (iind[0] = 0; iind[0] < A2.size(); ++iind[0])
    for (iind[1] = 0; iind[1] <= iind[0]; ++iind[1]) {
      A2.put(iind, 6 * count);
      cout << "Should be initialized to " << 6 * count << ":  result("
           << iind[0] << "," << iind[1] << ") = " << A2.get(iind) << endl;
      count++;
    }

  A2.resize(6);
  A2 = 25.7;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2 should be set to 25.7:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 += 6.6;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2 should be set to 32.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 -= 10;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2 should be set to 22.3:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 *= 10;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2 should be set to 223.0:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }
  A2 /= 100;
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2 should be set to 2.230:  result(" << j << "," << k
           << ") = " << A2(j, k) << endl;
    }

  A2.resize(3);
  A1.resize(3);
  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      A2(j, k) = 3.4 * (j + k);
      A1(j, k) = -2.7 * (j - 3 * k);
    }

  for (int j = 0; j < int(A2.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A2(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A2(j, k) << endl;
      cout << "A1(" << j << "," << k << ") should be " << -2.7 * (j - 3 * k)
           << ": result = " << A1(j, k) << endl;
    }

  LTMatrix A3;
  A3 = A2;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A3(" << j << "," << k << ") should be " << 3.4 * (j + k)
           << ": result = " << A3(j, k) << endl;
    }
  A3 += A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) - 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }
  A3 = A2;
  A3 -= A1;
  for (int j = 0; j < int(A3.size()); ++j)
    for (int k = 0; k <= j; ++k) {
      cout << "A3(" << j << "," << k << ") should be "
           << 3.4 * (j + k) + 2.7 * (j - 3 * k) << ": result = " << A3(j, k)
           << endl;
    }

  A2.resize(6);
  try {
    cout << "Trying to add mismatched size: should throw exception" << endl;
    A1 += A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }
  try {
    cout << "Trying to subtract mismatched size: should throw exception"
         << endl;
    A1 -= A2;
    cout << "ERROR" << endl;
  } catch (const std::exception& xp) {
    cout << "Exception caught" << endl;
  }

  cout << "simple tests done" << endl << endl;
}

template <typename T>
bool dotestLowerTriangularMatrix(const vector<T>& datacmp,
                                 const vector<unsigned int>& idata) {
  unsigned int successcount = 0;
  unsigned int failcount = 0;

  int firstsize = 6;
  LowerTriangularMatrix<T> H0(firstsize);
  for (int col = 0; col < firstsize; col++)
    for (int row = col; row < firstsize; row++)
      H0(row, col) = 100 * col + row;

  for (int row = firstsize - 1; row >= 0; row--)
    for (int col = 0; col <= row; col++)
      if (equal(T(100 * col + row), H0(row, col)))
        successcount++;
      else
        failcount++;

  cout << "cf=" << datacmp[0] << endl;
  LowerTriangularMatrix<T> H1(firstsize, datacmp[0]);
  for (int row = firstsize - 1; row >= 0; row--)
    for (int col = 0; col <= row; col++)
      if (equal(H1(row, col), datacmp[0]))
        successcount++;
      else
        failcount++;

  int secondsize = idata[3] * 128;
  cout << "secondsize = " << secondsize << endl;
  LowerTriangularMatrix<T> H2(secondsize);
  unsigned int count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      H2(row, col) = datacmp[count++];
      if (count > datacmp.size())
        count = 0;
    }

  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H2(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      H2.put(row, col, datacmp[count++]);
      if (count > datacmp.size())
        count = 0;
    }

  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H2.get(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  LowerTriangularMatrix<T> H3;
  H3 = H2;
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H3(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  LowerTriangularMatrix<T> H4(H2);
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H3(row, col), datacmp[count++]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H1;
  for (int row = firstsize - 1; row >= 0; row--)
    for (int col = 0; col <= row; col++)
      if (equal(H4(row, col), datacmp[0]))
        successcount++;
      else
        failcount++;

  H4.resize(4);
  H4 = datacmp[5];
  for (int col = 0; col < 4; col++)
    for (int row = col; row < 4; row++) {
      if (equal(H4(row, col), datacmp[5]))
        successcount++;
      else
        failcount++;
    }

  H4.clear();
  if (H4.size() == 0)
    successcount++;
  else
    failcount++;

  H3 += datacmp[8];
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H3(row, col), datacmp[count++] + datacmp[8]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 -= datacmp[6];
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H3(row, col), datacmp[count++] + datacmp[8] - datacmp[6]))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H3 *= datacmp[2];
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H3(row, col),
                datacmp[2] * (datacmp[count++] + datacmp[8] - datacmp[6])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  H4 = H3;
  H4 += H3;
  count = 0;
  for (int col = 0; col < secondsize; col++)
    for (int row = col; row < secondsize; row++) {
      if (equal(H4(row, col),
                2 * datacmp[2] * (datacmp[count++] + datacmp[8] - datacmp[6])))
        successcount++;
      else
        failcount++;
      if (count > datacmp.size())
        count = 0;
    }

  cout << endl;
  cout << " Number of successful tests = " << successcount << endl;
  cout << "     Number of FAILED tests = " << failcount << endl;
  return (failcount == 0);
}

void testLowerTriangularMatrix() {
  cout << "Starting LowerTriangularMatrix tests" << endl;

  dotestLTMatrixSimple();

  bool success = true;
  srand(time(NULL));
  unsigned int nintegers = 2056;
  vector<unsigned int> idata(nintegers);
  for (unsigned int k = 0; k < nintegers; ++k)
    idata[k] = rand() % 8 + 1;

  unsigned int ndata = 5048;
  vector<float> fdata(ndata);
  for (unsigned int k = 0; k < ndata; ++k)
    fdata[k] =
        1.434 * ((rand() % 4096 + 1) - 2047) - 0.54511 * ((rand() % 128) - 63);

  success = dotestLowerTriangularMatrix(fdata, idata);

  if (success)
    cout << "ALL LowerTriangularMatrix TESTS PASSED!!" << endl;
  else
    cout << "Some LowerTriangularMatrix tests FAILED" << endl;
}

void diagonalMatrix(RealSymmetricMatrix& mat, const double* diag, int n) {
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat(k, k) = diag[k];
}

void diagonalMatrix(ComplexHermitianMatrix& mat, const double* diag, int n) {
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat.put(k, k, complex<double>(diag[k], 0.0));
}

void diagonalMatrix(ComplexHermitianMatrix& mat, const RVector& lambda) {
  int n = lambda.size();
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat.put(k, k, complex<double>(lambda[k], 0.0));
}

void diagonalMatrix(RealSymmetricMatrix& mat, const RVector& lambda) {
  int n = lambda.size();
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat(k, k) = lambda[k];
}

void identity_matrix(ComplexHermitianMatrix& mat, int n) {
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat.put(k, k, complex<double>(1.0, 0.0));
}

void identity_matrix(RealSymmetricMatrix& mat, int n) {
  mat.resize(n);
  mat = 0.0;
  for (int k = 0; k < n; ++k)
    mat(k, k) = 1.0;
}

void herm_adj(CMatrix& out, const CMatrix& in) {
  out.resize(in.size(1), in.size(0));
  for (int row = 0; row < int(in.size(0)); row++)
    for (int col = 0; col < int(in.size(1)); col++)
      out(col, row) = conjugate(in(row, col));
}

void herm_adj(RMatrix& out, const RMatrix& in) {
  out.resize(in.size(1), in.size(0));
  for (int row = 0; row < int(in.size(0)); row++)
    for (int col = 0; col < int(in.size(1)); col++)
      out(col, row) = in(row, col);
}

void multiply(CMatrix& out, const CMatrix& inA, const CMatrix& inB) {
  if (inA.size(1) != inB.size(0))
    throw(std::invalid_argument("bad matrix multiplying sizes"));
  out.resize(inA.size(0), inB.size(1));
  for (int row = 0; row < int(out.size(0)); row++)
    for (int col = 0; col < int(out.size(1)); col++) {
      complex<double> z(0.0, 0.0);
      for (int k = 0; k < int(inA.size(1)); k++)
        z += inA(row, k) * inB(k, col);
      out(row, col) = z;
    }
}

void multiply(RMatrix& out, const RMatrix& inA, const RMatrix& inB) {
  if (inA.size(1) != inB.size(0))
    throw(std::invalid_argument("bad matrix multiplying sizes"));
  out.resize(inA.size(0), inB.size(1));
  for (int row = 0; row < int(out.size(0)); row++)
    for (int col = 0; col < int(out.size(1)); col++) {
      double z = 0.0;
      for (int k = 0; k < int(inA.size(1)); k++)
        z += inA(row, k) * inB(k, col);
      out(row, col) = z;
    }
}

CMatrix matrix_convert(const ComplexHermitianMatrix& mat) {
  CMatrix temp(mat.size(), mat.size());
  for (int i = 0; i < int(mat.size()); i++)
    for (int j = 0; j < int(mat.size()); j++)
      temp(i, j) = mat(i, j);
  return temp;
}

RMatrix matrix_convert(const RealSymmetricMatrix& mat) {
  RMatrix temp(mat.size(), mat.size());
  for (int i = 0; i < int(mat.size()); i++)
    for (int j = 0; j < int(mat.size()); j++)
      temp(i, j) = mat(i, j);
  return temp;
}

bool is_equal(const CMatrix& mat1, const CMatrix& mat2) {
  if (mat1.size(0) != mat2.size(0))
    return false;
  if (mat1.size(1) != mat2.size(1))
    return false;
  for (int k = 0; k < int(mat1.size(0)); k++)
    for (int l = 0; l < int(mat1.size(1)); l++)
      if (std::abs(mat1(k, l) - mat2(k, l)) > 1e-10)
        return false;
  return true;
}

bool is_equal(const RMatrix& mat1, const RMatrix& mat2) {
  if (mat1.size(0) != mat2.size(0))
    return false;
  if (mat1.size(1) != mat2.size(1))
    return false;
  for (int k = 0; k < int(mat1.size(0)); k++)
    for (int l = 0; l < int(mat1.size(1)); l++)
      if (std::abs(mat1(k, l) - mat2(k, l)) > 1e-10)
        return false;
  return true;
}

void mat_print(const CMatrix& mat, const string& matname) {
  cout << endl;
  for (int i = 0; i < int(mat.size(0)); i++)
    for (int j = 0; j < int(mat.size(1)); j++)
      cout << matname << "(" << i << "," << j << ") = " << mat(i, j) << endl;
  cout << endl;
}

void mat_print(const RMatrix& mat, const string& matname) {
  cout << endl;
  for (int i = 0; i < int(mat.size(0)); i++)
    for (int j = 0; j < int(mat.size(1)); j++)
      cout << matname << "(" << i << "," << j << ") = " << mat(i, j) << endl;
  cout << endl;
}

void testMatrix(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMatrix") == 0)
    return;

  cout << "TestMatrix" << endl;
  cout << endl
       << endl
       << "*******testVector**************************************************"
       << endl
       << endl;
  testVector();
  cout << endl
       << endl
       << "*******testMatrix**************************************************"
       << endl
       << endl;
  testMatrix();
  cout << endl
       << endl
       << "*******testLowerTriangularMatrix************************************"
          "**************"
       << endl
       << endl;
  testLowerTriangularMatrix();
  cout << endl
       << endl
       << "*******testRealSymmetricMatrix**************************************"
          "************"
       << endl
       << endl;
  testRealSymmetricMatrix();
  cout << endl
       << endl
       << "*******testComplexHermitianMatrix***********************************"
          "***************"
       << endl
       << endl;
  testComplexHermitianMatrix();
  cout << endl
       << endl
       << "*******testRealDiagonalize******************************************"
          "********"
       << endl
       << endl;
}

// ******************************************************************************
