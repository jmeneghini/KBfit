#ifndef TASK_UTILS_H
#define TASK_UTILS_H

#include "matrix.h"

// ***************************************************************************************

//  Takes a Hermitian matrix "H" and returns the eigenvalues in
//  ascending order in "eigvals" and the associated eigenvectors
//  in the columns of "eigvecs".  Throws an exception if fails.
//  Versions for only the eigenvalues are also available.

class Diagonalizer {
public:
  Diagonalizer() {}
  ~Diagonalizer() {}

  void getEigenvalues(const RealSymmetricMatrix& H, RVector& eigvals);
  void getEigenvectors(const RealSymmetricMatrix& H, RVector& eigvals,
                       RMatrix& eigvecs);
  void getEigenvalues(const ComplexHermitianMatrix& H, RVector& eigvals);
  void getEigenvectors(const ComplexHermitianMatrix& H, RVector& eigvals,
                       CMatrix& eigvecs);

private:
  void diagonalize(const RealSymmetricMatrix& H, RVector& eigvals,
                   RMatrix& eigvecs, bool calceigvecs);
  void diagonalize(const ComplexHermitianMatrix& H, RVector& eigvals,
                   CMatrix& eigvecs, bool calceigvecs);
};

// ****************************************************************

//  Given a positive definite real symmetric matrix "A",
//  an object of this class can construct its Cholesky decomposition
//                 A = L * transpose(L)
//  where L is lower triangular, or the Cholesky decomposition of the
//  inverse of A:
//                 A^(-1) = transpose(L) * L
//  where L is lower triangular.  Throws an exception if not
//  successful.

class CholeskyDecomposer {
public:
  CholeskyDecomposer() {}
  ~CholeskyDecomposer() {}

  void getCholesky(const RealSymmetricMatrix& A,
                   LowerTriangularMatrix<double>& L);
  void getCholeskyOfInverse(const RealSymmetricMatrix& A,
                            LowerTriangularMatrix<double>& L);
};

// **************************************************************
#endif
