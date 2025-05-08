#ifndef UTILS_H
#define UTILS_H

#include "matrix.h"

typedef std::vector<double> Rvector;
typedef std::vector<std::complex<double>> Cvector;

// ***************************************************************************************

//  Takes a Hermitian matrix "H" and computes its inverse in "Hinv".

void calcMatrixInverse(const ComplexHermitianMatrix& H,
                       ComplexHermitianMatrix& Hinv);

//  Takes a general complex matrix "M" and computes its inverse in "Minv".

void calcMatrixInverse(const CMatrix& M, CMatrix& Minv);

//  Computes the Cayley transformation of a Hermitian matrix
//        computes  -d*(d+i*H)*(-d+i*H)^(-1)
//  with   d = (Mform) ? 1.0 : -1.0
//     Mform true is the form for H = K matrix,
//     Mform false is the form for H = K inverse

void calcCayleyTransformMatrix(const ComplexHermitianMatrix& H,
                               CMatrix& CTMatrix, bool Mform);

void calcCayleyTransformMatrix(const RealSymmetricMatrix& H, CMatrix& CTMatrix,
                               bool Mform);

// ***************************************************************************************

//  Takes a Hermitian matrix "H" and returns the eigenvalues in
//  ascending order in "eigvals" and the associated eigenvectors
//  in the columns of "eigvecs".  Throws an exception if fails.
//  Versions for only the eigenvalues are also available.
//  The complex eigenvalues of a complex non-Hermitian matrix
//  can also be found.

class Diagonalizer {
public:
  Diagonalizer() {}
  ~Diagonalizer() {}

  void getEigenvalues(const RealSymmetricMatrix& H, Rvector& eigvals);
  void getEigenvectors(const RealSymmetricMatrix& H, Rvector& eigvals,
                       RMatrix& eigvecs);
  void getEigenvalues(const ComplexHermitianMatrix& H, Rvector& eigvals);
  void getEigenvectors(const ComplexHermitianMatrix& H, Rvector& eigvals,
                       CMatrix& eigvecs);
  void getEigenvalues(const CMatrix& M, Cvector& eigvals);

private:
  void diagonalize(const RealSymmetricMatrix& H, Rvector& eigvals,
                   RMatrix& eigvecs, bool calceigvecs);
  void diagonalize(const ComplexHermitianMatrix& H, Rvector& eigvals,
                   CMatrix& eigvecs, bool calceigvecs);
};

// ****************************************************************

//  Evaluates
//     det(Q)  [returned as a product of complex numbers]
//     Omega(mu,Q) = det(Q) / det(sqrt(mu^2+Q*Qdagger))
//  The determinant is evaluated using LU decomposition, so
//  the values returned are NOT the eigenvalues.  The Omega
//  function is evaluated using LU decomposition in the
//  numerator, and Cholesky decomposition in the denominator.

class DeterminantCalculator {
public:
  DeterminantCalculator() {}
  ~DeterminantCalculator() {}

  void getDeterminantAsProduct(const CMatrix& M, Cvector& DetProd);
  std::complex<double> getOmega(double mu, const CMatrix& M);

private:
  std::complex<double> get_omega(std::vector<double>& mat,
                                 std::vector<double>& bmatf, int n);
};

// ****************************************************************
#endif
