#include "task_utils.h"
using namespace std;

// ***************************************************************************************

// Prototypes of routines in LAPACK library--to call Fortran
// routines from a C++ program, use extern "C" to tell the
// compiler that the external routine is a C routine; then
// add an underscore to the end of the routine name since
// the routine is in Fortran.  All parameters must be passed
// as pointers and two-dimensional arrays must follow
// Fortran conventions.

extern "C" {
void dpotrf_(char*, int*, double*, int*, int*);
void dpotri_(char*, int*, double*, int*, int*);
void dsygv_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda,
            double* b, int* ldb, double* w, double* work, int* lwork,
            int* info);
void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, int* info);
void zheev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, double* rwork, int* info);
}

// ***************************************************************

//  Takes a Hermitian matrix "H" and returns the eigenvalues in
//  ascending order in "eigvals" and the associated eigenvectors
//  in the columns of "eigvecs".  Throws an exception if fails.

void Diagonalizer::diagonalize(const RealSymmetricMatrix& H, RVector& eigvals,
                               RMatrix& eigvecs, bool calceigvecs) {
  int n = H.size();
  if (n == 0) {
    eigvals.clear();
    eigvecs.clear();
    return;
  }
  int lwork = 5 * n;
  RVector work(lwork);
  eigvals.resize(n);
  int info;
  char jobz = (calceigvecs) ? 'V' : 'N';
  char uplo = 'U';

  // load H (upper triangle) into matf fortran format
  //    (column major; row index changes fastest)
  vector<double> matf(n * n);
  for (int col = 0; col < n; ++col)
    for (int row = 0; row <= col; ++row)
      matf[row + n * col] = H(row, col);

  dsyev_(&jobz, &uplo, &n, &matf[0], &n, &eigvals[0], &work[0], &lwork, &info);
  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in diagonalize"));
  } else if (info > 0) {
    throw(std::invalid_argument(" no convergence in diagonalize"));
  }

  if (calceigvecs) {
    eigvecs.resize(n, n);
    for (int col = 0; col < n; ++col)
      for (int row = 0; row < n; ++row)
        eigvecs(row, col) = matf[row + n * col];
  }
}

void Diagonalizer::getEigenvectors(const RealSymmetricMatrix& H,
                                   RVector& eigvals, RMatrix& eigvecs) {
  diagonalize(H, eigvals, eigvecs, true);
}

void Diagonalizer::getEigenvalues(const RealSymmetricMatrix& H,
                                  RVector& eigvals) {
  RMatrix eigvecs;
  diagonalize(H, eigvals, eigvecs, false);
}

void Diagonalizer::diagonalize(const ComplexHermitianMatrix& H,
                               RVector& eigvals, CMatrix& eigvecs,
                               bool calceigvecs) {
  int n = H.size();
  if (n == 0) {
    eigvals.clear();
    eigvecs.clear();
    return;
  }
  int lwork = 4 * n;
  RVector work(2 * lwork);
  RVector rwork(3 * n);
  eigvals.resize(n);
  int info;
  char jobz = (calceigvecs) ? 'V' : 'N';
  char uplo = 'U';

  // load H (upper triangle) into matf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> matf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row <= col; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = H(row, col);
      matf[index] = z.real();
      matf[index + 1] = z.imag();
    }

  zheev_(&jobz, &uplo, &n, &matf[0], &n, &eigvals[0], &work[0], &lwork,
         &rwork[0], &info);
  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in diagonalize"));
  } else if (info > 0) {
    throw(std::invalid_argument(" no convergence in diagonalize"));
  }

  // cout << "optimal lwork = "<<work[0]<<endl;

  if (calceigvecs) {
    eigvecs.resize(n, n);
    for (int col = 0; col < n; col++)
      for (int row = 0; row < n; row++) {
        int index = 2 * (row + n * col);
        eigvecs(row, col) = complex<double>(matf[index], matf[index + 1]);
      }
  }
}

void Diagonalizer::getEigenvectors(const ComplexHermitianMatrix& H,
                                   RVector& eigvals, CMatrix& eigvecs) {
  diagonalize(H, eigvals, eigvecs, true);
}

void Diagonalizer::getEigenvalues(const ComplexHermitianMatrix& H,
                                  RVector& eigvals) {
  CMatrix eigvecs;
  diagonalize(H, eigvals, eigvecs, false);
}

// ****************************************************************

// Given a positive definite real symmetric matrix "A",
// this routine constructs its Cholesky decomposition
//                 A = L * transpose(L)
// where L is lower triangular.  Throws an exception if not successful.

void CholeskyDecomposer::getCholesky(const RealSymmetricMatrix& A,
                                     LowerTriangularMatrix<double>& L) {
  int n = A.size();
  if (n == 0) {
    L.clear();
    return;
  }
  int info;
  char uplo = 'L';

  // load A into lower triangle of mata in fortran format
  //    (column major; row index changes fastest)
  vector<double> mata(n * n);
  for (int row = 0; row < n; ++row)
    for (int col = 0; col <= row; ++col)
      mata[row + n * col] = A(row, col);

  dpotrf_(&uplo, &n, &mata[0], &n, &info);
  if (info < 0) {
    L.clear();
    throw(std::invalid_argument(" bad arguments in cholesky"));
  } else if (info > 0) {
    L.clear();
    throw(std::invalid_argument(" matrix not positive definite in cholesky"));
  }

  L.resize(n);
  for (int row = 0; row < n; ++row)
    for (int col = 0; col <= row; ++col)
      L(row, col) = mata[row + n * col];
}

// *************************************************************

// Given a positive definite real symmetric matrix "A",
// this routine constructs the Cholesky decomposition of the
// inverse of A:
//                 A^(-1) = transpose(L) * L
// where L is lower triangular. Throws an exception if not successful.

void CholeskyDecomposer::getCholeskyOfInverse(
    const RealSymmetricMatrix& A, LowerTriangularMatrix<double>& L) {
  try {
    getCholesky(A, L);
  } catch (const std::exception& errmsg) {
    throw(std::invalid_argument(string("Failure in cholesky_of_inverse: ") +
                                string(errmsg.what())));
  }
  int n = A.size();
  for (int i = 0; i < n; i++) {
    L(i, i) = 1.0 / L(i, i);
    for (int j = i + 1; j < n; j++) {
      double sum = 0.0;
      for (int k = i; k < j; k++)
        sum -= L(j, k) * L(k, i);
      L(j, i) = sum / L(j, j);
    }
  }
}

// ************************************************************
