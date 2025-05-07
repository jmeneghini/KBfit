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

ComplexHermitianMatrix& ComplexHermitianMatrix::invertRoot() {
  int n = static_cast<int>(m_size);

  RVector eigvals(n);
  RVector inv_root_eigvals(n);
  CMatrix eigvecs;
  Diagonalizer D;

  // Compute the eigenvalues and eigenvectors of the matrix.
  D.getEigenvectors(*this, eigvals, eigvecs);

  //  // Print the eigenvecs
  //  std::cout << "Eigenvecs: " << std::endl;
  //  for (int i = 0; i < n; ++i)
  //  {
  //     for (int j = 0; j < n; ++j)
  //     {
  //        std::cout << eigvecs(i, j) << " ";
  //     }
  //     std::cout << std::endl;
  //  }

  // Check if eigenvalues are positive. If so, replace with inverse square root.
  for (int i = 0; i < n; ++i) {
    inv_root_eigvals[i] = 1.0 / eigvals[i];
  }

  // Reconstruct A^(-1/2) = U diag(1/sqrt(eigvals)) U^dag -> A^(-1/2)_ij = sum_k
  // U_ik (1/sqrt(eigvals_k)) U_jk^*
  std::vector<std::complex<double>> result(0.5 * n * (n + 1));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      for (int k = 0; k < n; ++k) {
        result[get_index(i, j)] =
            eigvecs(i, k) * conjugate(eigvecs(j, k)) * inv_root_eigvals[k];
      }
    }
  }

  m_store = result;

  //  std::cout << "InvertRoot: " << std::endl;
  //  for (int i = 0; i < n; ++i)
  //  {
  //     for (int j = 0; j < n; ++j)
  //     {
  //        std::cout << m_store[get_index(i, j)] << " ";
  //     }
  //     std::cout << std::endl;
  //  }

  return *this;
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
