#ifndef HELPER_H
#define HELPER_H
#include "matrix.h"

namespace TestHelper {
  inline bool approx_equal_complex(const std::complex<double>& a, const std::complex<double>& b, double epsilon = 1e-6);
  inline bool approx_equal_cmatrix(const CMatrix& A, const CMatrix& B, double epsilon = 1e-6);
  inline bool is_unitary(const CMatrix& M, double epsilon = 1e-6);
}

bool TestHelper::approx_equal_complex(const std::complex<double>& a, const std::complex<double>& b, double epsilon) {
  return std::abs(a.real() - b.real()) < epsilon && std::abs(a.imag() - b.imag()) < epsilon;
}

// Helper for comparing CMatrix objects
bool TestHelper::approx_equal_cmatrix(const CMatrix& A, const CMatrix& B, double epsilon) {
  if (A.size(0) != B.size(0) || A.size(1) != B.size(1)) {
    // Using doctest::Context().cout() for output within tests if needed for debugging
    // doctest::Context().cout() << "Matrix size mismatch: A(" << A.size(0) << "," << A.size(1)
    //                           << ") vs B(" << B.size(0) << "," << B.size(1) << ")" << std::endl;
    return false;
  }
  for (unsigned r = 0; r < A.size(0); ++r) {
    for (unsigned c = 0; c < A.size(1); ++c) {
      // Assuming CMatrix has an operator() or a get() method for element access
      if (!approx_equal_complex(A(r, c), B(r, c), epsilon)) {
        // doctest::Context().cout() << "Mismatch at (" << r << "," << c << "): A=" << A(r,c) << ", B=" << B(r,c) << std::endl;
        return false;
      }
    }
  }
  return true;
}

bool TestHelper::is_unitary(const CMatrix& M, double epsilon) {
  if (M.size(0) != M.size(1) || M.size(0) == 0) { // Must be square and non-empty
    // doctest::Context().cout() << "Matrix is not square or is empty for unitary check." << std::endl;
    return false;
  }
  unsigned n = M.size(0);
  CMatrix prod(n, n); // To store M * M_dagger

  // Calculate M * M_dagger
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      std::complex<double> sum_val(0.0, 0.0);
      for (unsigned k = 0; k < n; ++k) {
        sum_val += M(i, k) * std::conj(M(j, k)); // (M * M_dagger)_ij = sum_k M_ik * conj(M_jk)
      }
      prod.put(i, j, sum_val);
    }
  }

  // Create an identity matrix of the same size
  CMatrix identity(n, n);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      if (i == j) {
        identity.put(i, j, std::complex<double>(1.0, 0.0));
      } else {
        identity.put(i, j, std::complex<double>(0.0, 0.0));
      }
    }
  }
  return approx_equal_cmatrix(prod, identity, epsilon);
}

#endif



