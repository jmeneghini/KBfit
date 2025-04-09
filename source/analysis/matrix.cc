#include "matrix.h"
#include "task_utils.h"
using namespace std;

// **************************************************************

double conjugate(const double& x) { return x; }

float conjugate(const float& x) { return x; }

std::complex<double> conjugate(const std::complex<double>& z) {
  return conj(z);
}

std::complex<float> conjugate(const std::complex<float>& z) { return conj(z); }

double realpart(const double& x) { return x; }

double realpart(const float& x) { return double(x); }

double imaginarypart(const double& x) { return 0.0; }

double imaginarypart(const float& x) { return 0.0; }

double realpart(const complex<double>& z) { return real(z); }

double realpart(const complex<float>& z) { return double(real(z)); }

double imaginarypart(const complex<double>& z) { return imag(z); }

double imaginarypart(const complex<float>& z) { return double(imag(z)); }

double sqr(const double& x) { return x * x; }

double sqr(const std::complex<double>& z) { return std::norm(z); }

// *********************************************************************

RealSymmetricMatrix::RealSymmetricMatrix() { m_size = 0; }

RealSymmetricMatrix::RealSymmetricMatrix(int insize)
    : m_store(set_size(insize)) {}

RealSymmetricMatrix::RealSymmetricMatrix(uint insize)
    : m_store(set_size(insize)) {}

RealSymmetricMatrix::RealSymmetricMatrix(int insize, double initial_value)
    : m_store(set_size(insize), initial_value) {}

RealSymmetricMatrix::RealSymmetricMatrix(uint insize, double initial_value)
    : m_store(set_size(insize), initial_value) {}

RealSymmetricMatrix::RealSymmetricMatrix(const RealSymmetricMatrix& incoming)
    : m_store(incoming.m_store), m_size(incoming.m_size) {}

RealSymmetricMatrix::~RealSymmetricMatrix() {}

RealSymmetricMatrix& RealSymmetricMatrix::clear() {
  m_store.clear();
  m_size = 0;
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] = val;
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator+=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] += val;
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator-=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] -= val;
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator*=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] *= val;
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::operator/=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] /= val;
  return *this;
}

RealSymmetricMatrix&
RealSymmetricMatrix::operator=(const RealSymmetricMatrix& incoming) {
  if (this == &incoming)
    return *this;
  m_store = incoming.m_store;
  m_size = incoming.m_size;
  return *this;
}

RealSymmetricMatrix&
RealSymmetricMatrix::operator+=(const RealSymmetricMatrix& incoming) {
#ifdef SAFETY_FLAG
  if (m_size != incoming.m_size)
    throw(std::invalid_argument("RealSymmetricMatrix size mismatch"));
#endif
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] += incoming.m_store[i];
  return *this;
}

RealSymmetricMatrix&
RealSymmetricMatrix::operator-=(const RealSymmetricMatrix& incoming) {
#ifdef SAFETY_FLAG
  if (m_size != incoming.m_size)
    throw(std::invalid_argument("RealSymmetricMatrix size mismatch"));
#endif
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] -= incoming.m_store[i];
  return *this;
}

void RealSymmetricMatrix::modifyEigenvalues(std::complex<double> (*func)(double), Matrix<std::complex<double>>& out_matrix) const {
  int n = static_cast<int>(m_size);

  RVector eigvals(n);
  std::vector<complex<double>> modified_eigvals(n);
  RMatrix eigvecs;
  Diagonalizer D;

  // Compute the eigenvalues and eigenvectors of the matrix.
  D.getEigenvectors(*this, eigvals, eigvecs);

  for (int i = 0; i < n; ++i) {
    modified_eigvals[i] = func(eigvals[i]);
  }


  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      complex<double> element(0.0, 0.0);
      for (int k = 0; k < n; ++k) {
        element += eigvecs(i, k) * eigvecs(j, k) * modified_eigvals[k];
      }
      out_matrix.put(i, j, element);
    }
  }
}

RealSymmetricMatrix& RealSymmetricMatrix::resize() { return clear(); }

RealSymmetricMatrix& RealSymmetricMatrix::resize(int insize) {
  m_store.resize(set_size(insize));
  return *this;
}

RealSymmetricMatrix& RealSymmetricMatrix::resize(uint insize) {
  m_store.resize(set_size(insize));
  return *this;
}

// **************************************************************

ComplexHermitianMatrix::ComplexHermitianMatrix() { m_size = 0; }

ComplexHermitianMatrix::ComplexHermitianMatrix(int insize)
    : m_store(set_size(insize)) {}

ComplexHermitianMatrix::ComplexHermitianMatrix(uint insize)
    : m_store(set_size(insize)) {}

ComplexHermitianMatrix::ComplexHermitianMatrix(int insize, double initial_value)
    : m_store(set_size(insize), std::complex<double>(initial_value, 0.0)) {}

ComplexHermitianMatrix::ComplexHermitianMatrix(uint insize,
                                               double initial_value)
    : m_store(set_size(insize), std::complex<double>(initial_value, 0.0)) {}

ComplexHermitianMatrix::ComplexHermitianMatrix(
    const ComplexHermitianMatrix& incoming)
    : m_store(incoming.m_store), m_size(incoming.m_size) {}

ComplexHermitianMatrix::~ComplexHermitianMatrix() {}

ComplexHermitianMatrix& ComplexHermitianMatrix::clear() {
  m_store.clear();
  m_size = 0;
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] = val;
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator+=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] += val;
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator-=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] -= val;
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator*=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] *= val;
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::operator/=(double val) {
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] /= val;
  return *this;
}

ComplexHermitianMatrix&
ComplexHermitianMatrix::operator=(const ComplexHermitianMatrix& incoming) {
  if (this == &incoming)
    return *this;
  m_store = incoming.m_store;
  m_size = incoming.m_size;
  return *this;
}

ComplexHermitianMatrix&
ComplexHermitianMatrix::operator+=(const ComplexHermitianMatrix& incoming) {
#ifdef SAFETY_FLAG
  if (m_size != incoming.m_size)
    throw(std::invalid_argument("ComplexHermitianMatrix size mismatch"));
#endif
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] += incoming.m_store[i];
  return *this;
}

ComplexHermitianMatrix&
ComplexHermitianMatrix::operator-=(const ComplexHermitianMatrix& incoming) {
#ifdef SAFETY_FLAG
  if (m_size != incoming.m_size)
    throw(std::invalid_argument("ComplexHermitianMatrix size mismatch"));
#endif
  for (uint i = 0; i < m_store.size(); ++i)
    m_store[i] -= incoming.m_store[i];
  return *this;
}

void ComplexHermitianMatrix::modifyEigenvalues(std::complex<double> (*func)(double), Matrix<std::complex<double>>& out_matrix) const {
  int n = static_cast<int>(m_size);

  RVector eigvals(n);
  std::vector<complex<double>> modified_eigvals(n);
  CMatrix eigvecs;
  Diagonalizer D;

  // Compute the eigenvalues and eigenvectors of the matrix.
  D.getEigenvectors(*this, eigvals, eigvecs);

  for (int i = 0; i < n; ++i) {
    modified_eigvals[i] = func(eigvals[i]);
  }


  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      complex<double> element(0.0, 0.0);
      for (int k = 0; k < n; ++k) {
        element += eigvecs(i, k) * conjugate(eigvecs(j, k)) * modified_eigvals[k];
      }
      out_matrix.put(i, j, element);
    }
  }
}

ComplexHermitianMatrix& ComplexHermitianMatrix::resize() { return clear(); }

ComplexHermitianMatrix& ComplexHermitianMatrix::resize(int insize) {
  m_store.resize(set_size(insize));
  return *this;
}

ComplexHermitianMatrix& ComplexHermitianMatrix::resize(uint insize) {
  m_store.resize(set_size(insize));
  return *this;
}

// ***************************************************************************

RVector copyRVectorWithoutZerothElement(const RVector& in) {
  uint n = in.size() - 1;
  RVector buffer(n);
  for (uint k = 0; k < n; k++)
    buffer[k] = in[k + 1];
  return buffer;
}

// **************************************************************
