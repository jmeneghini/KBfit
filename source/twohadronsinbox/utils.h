#ifndef UTILS_H
#define UTILS_H

#include "matrix.h"

typedef std::vector<double> Rvector;
typedef std::vector<std::complex<double>> Cvector;

// ***************************************************************************************

//  Takes a Hermitian matrix "H" and computes its inverse in "Hinv".

void calcMatrixInverse(const ComplexHermitianMatrix& H,
                       ComplexHermitianMatrix& Hinv);

//  Calculates the phase angle of the complex number "z"
//  in the range [0, 2*PI) and returns it in "phase".
//  This is in contrast to the standard atan2 function which
//  returns the angle in the range (-PI, PI].
double getPhaseAngle(const std::complex<double>& z);

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
  void getEigenvectors(const CMatrix& M, Cvector& eigvals, CMatrix& eigvecs);

private:
  void diagonalize(const RealSymmetricMatrix& H, Rvector& eigvals,
                   RMatrix& eigvecs, bool calceigvecs);
  void diagonalize(const ComplexHermitianMatrix& H, Rvector& eigvals,
                   CMatrix& eigvecs, bool calceigvecs);
};

// ****************************************************************

//  Evaluates
//              det(Q)^(1/Ndet)  and
//              Omega(mu,Q) = det(Q) / det(sqrt(mu^2+Q*Qdagger))
//  where Ndet is a positive, odd integer, and
//
//    (a) either Q = M where M is Hermitian,
//    (b) or  Q = 1 - A*B  where A real symmetric
//          and B is Hermitian.
//  These quantities must be real.  The Ndet root avoids overflow
//  and underflow.

class RealDeterminantRoot {
public:
  RealDeterminantRoot() {}
  ~RealDeterminantRoot() {}

  double getDeterminantOddRoot(const ComplexHermitianMatrix& M, uint Ndet);

  double getDeterminantOddRoot(const RealSymmetricMatrix& A,
                               const ComplexHermitianMatrix& B, uint Ndet);

  double getOmega(double mu, const RVector& eigenvalues);

  double getOmega(double mu, const CVector& eigenvalues, double& imag_part);

  double getOmega(double mu, const ComplexHermitianMatrix& M);

  double getOmega(double mu, const RealSymmetricMatrix& A,
                  const ComplexHermitianMatrix& B);

private:
  double get_det_odd_root(std::vector<double>& matfortran, int n, uint Ndet);

  double get_omega(std::vector<double>& mat, std::vector<double>& bmatf, int n);
};

// ****************************************************************
#endif
