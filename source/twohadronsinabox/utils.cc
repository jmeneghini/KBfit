#include "utils.h"
using namespace std;

// *************************************************************************

// Prototypes of routines in LAPACK library--to call Fortran
// routines from a C++ program, use extern "C" to tell the
// compiler that the external routine is a C routine; then
// add an underscore to the end of the routine name since
// the routine is in Fortran.  All parameters must be passed
// as pointers and two-dimensional arrays must follow
// Fortran conventions.

extern "C" {
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, int* info);
void zheev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
            double* work, int* lwork, double* rwork, int* info);
void zgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, double* w,
            double* vl, int* ldvl, double* vr, int* ldvr, double* work,
            int* lwork, double* rwork, int* info);
void zgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void zpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
void zgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork,
             int* info);
}

// ***************************************************************

//  Takes a Hermitian matrix "H" and computes its inverse in "Hinv".

void calcMatrixInverse(const ComplexHermitianMatrix& H,
                       ComplexHermitianMatrix& Hinv) {
  int n = H.size();
  Hinv.clear();
  Hinv.resize(n);
  if (n == 0)
    return;
  if (n == 1) {
    Hinv.put(0, 0, 1.0 / H(0, 0));
    return;
  }
  // load H (entire matrix) into matf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> matf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = H(row, col);
      matf[index] = z.real();
      matf[index + 1] = z.imag();
    }
  vector<int> ipiv(n);
  int info;
  zgetrf_(&n, &n, &matf[0], &n, &ipiv[0], &info);
  if (info < 0) {
    throw(std::invalid_argument(" calcMatrixInverse failed"));
  }
  vector<double> work(2 * n);

  zgetri_(&n, &matf[0], &n, &ipiv[0], &work[0], &n, &info);
  if (info < 0) {
    throw(std::invalid_argument(" calcMatrixInverse failed"));
  }

  int index = 0;
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row <= col; ++row) {
      if (row != col)
        Hinv.put(row, col, complex<double>(matf[index], matf[index + 1]));
      else
        Hinv.put(row, col, complex<double>(matf[index], 0.0));
      index += 2;
    }
    for (int row = col + 1; row < n; ++row)
      index += 2;
  }
}

// ***************************************************************

//  Takes a general square complex matrix "H" and
//  computes its inverse in "Hinv".

void calcMatrixInverse(const CMatrix& H, CMatrix& Hinv) {
  int n = H.size(0);
  if (H.size(1) != uint(n)) {
    throw(std::invalid_argument("calcMatrixInverse requires a square matrix"));
  }
  Hinv.clear();
  Hinv.resize(n, n);
  if (n == 0)
    return;
  if (n == 1) {
    Hinv.put(0, 0, 1.0 / H(0, 0));
    return;
  }
  // load H (entire matrix) into matf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> matf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = H(row, col);
      matf[index] = z.real();
      matf[index + 1] = z.imag();
    }
  vector<int> ipiv(n);
  int info;
  zgetrf_(&n, &n, &matf[0], &n, &ipiv[0], &info);
  if (info < 0) {
    throw(std::invalid_argument(" calcMatrixInverse failed"));
  }
  vector<double> work(2 * n);

  zgetri_(&n, &matf[0], &n, &ipiv[0], &work[0], &n, &info);
  if (info < 0) {
    throw(std::invalid_argument(" calcMatrixInverse failed"));
  }

  int index = 0;
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      Hinv.put(row, col, complex<double>(matf[index], matf[index + 1]));
      index += 2;
    }
  }
}

// ***************************************************************

//  Computes the Cayley transformation of a Hermitian matrix
//        computes  -d*(d+i*H)*(-d+i*H)^(-1)
//  with   d = (Mform) ? 1.0 : -1.0
//     Mform true is the form for H = K matrix,
//     Mform false is the form for H = K inverse

void calcCayleyTransformMatrix(const ComplexHermitianMatrix& H,
                               CMatrix& CTMatrix, bool Mform) {
  CTMatrix.clear();
  int n = H.size();
  if (n == 0)
    return;
  double d = (Mform) ? 1.0 : -1.0;
  CMatrix A(n, n);
  complex<double> z;
  // form A = -d+i*M
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < col; ++row) {
      z = H(row, col);
      A.put(row, col, complex<double>(-z.imag(), z.real()));
      A.put(col, row, complex<double>(z.imag(), z.real()));
    }
    z = H(col, col);
    A.put(col, col, complex<double>(-d - z.imag(), z.real()));
  }
  CMatrix B;
  // compute B = (-1+i*M)^(-1)
  calcMatrixInverse(A, B);
  // form A = 1 + i*M = 2+A
  double dd = 2.0 * d;
  for (int col = 0; col < n; ++col) {
    z = A(col, col);
    A.put(col, col, complex<double>(dd + z.real(), z.imag()));
  }
  // multiply A*B
  CTMatrix.resize(n, n);
  dd = -d;
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      z = complex<double>(0.0, 0.0);
      for (int ind = 0; ind < n; ++ind) {
        z += A(row, ind) * B(ind, col);
      }
      CTMatrix.put(row, col, dd * z);
    }
  }
}

void calcCayleyTransformMatrix(const RealSymmetricMatrix& H, CMatrix& CTMatrix,
                               bool Mform) {
  CTMatrix.clear();
  int n = H.size();
  if (n == 0)
    return;
  double d = (Mform) ? 1.0 : -1.0;
  CMatrix A(n, n);
  complex<double> z;
  double r;
  // form A = -d+i*M
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < col; ++row) {
      r = H(row, col);
      A.put(row, col, complex<double>(0.0, r));
      A.put(col, row, complex<double>(0.0, r));
    }
    r = H(col, col);
    A.put(col, col, complex<double>(-d, r));
  }
  CMatrix B;
  // compute B = (-1+i*M)^(-1)
  calcMatrixInverse(A, B);
  // form A = 1 + i*M = 2+A
  double dd = 2.0 * d;
  for (int col = 0; col < n; ++col) {
    z = A(col, col);
    A.put(col, col, complex<double>(dd + z.real(), z.imag()));
  }
  // multiply A*B
  CTMatrix.resize(n, n);
  dd = -d;
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      z = complex<double>(0.0, 0.0);
      for (int ind = 0; ind < n; ++ind) {
        z += A(row, ind) * B(ind, col);
      }
      CTMatrix.put(row, col, dd * z);
    }
  }
}

// ***************************************************************

//  Takes a Hermitian matrix "H" and returns the eigenvalues in
//  ascending order in "eigvals" and the associated eigenvectors
//  in the columns of "eigvecs".  Throws an exception if fails.

void Diagonalizer::diagonalize(const RealSymmetricMatrix& H, Rvector& eigvals,
                               RMatrix& eigvecs, bool calceigvecs) {
  int n = H.size();
  if (n == 0) {
    eigvals.clear();
    eigvecs.clear();
    return;
  }
  int lwork = 5 * n;
  Rvector work(lwork);
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
                                   Rvector& eigvals, RMatrix& eigvecs) {
  diagonalize(H, eigvals, eigvecs, true);
}

void Diagonalizer::getEigenvalues(const RealSymmetricMatrix& H,
                                  Rvector& eigvals) {
  RMatrix eigvecs;
  diagonalize(H, eigvals, eigvecs, false);
}

void Diagonalizer::diagonalize(const ComplexHermitianMatrix& H,
                               Rvector& eigvals, CMatrix& eigvecs,
                               bool calceigvecs) {
  int n = H.size();
  if (n == 0) {
    eigvals.clear();
    eigvecs.clear();
    return;
  }
  int lwork = 4 * n;
  Rvector work(2 * lwork);
  Rvector rwork(3 * n);
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
                                   Rvector& eigvals, CMatrix& eigvecs) {
  diagonalize(H, eigvals, eigvecs, true);
}

void Diagonalizer::getEigenvalues(const ComplexHermitianMatrix& H,
                                  Rvector& eigvals) {
  CMatrix eigvecs;
  diagonalize(H, eigvals, eigvecs, false);
}

void Diagonalizer::getEigenvalues(const CMatrix& M, Cvector& eigvals) {
  int n = M.size(0);
  if (int(M.size(1)) != n)
    throw(std::invalid_argument("Must be square matrix to get eigenvalues"));
  if (n == 0) {
    eigvals.clear();
    return;
  }

  int lwork = 4 * n;
  Rvector work(2 * lwork);
  Rvector rwork(2 * n);
  Rvector lambda(2 * n);
  int info;
  char jobvl = 'N';
  char jobvr = 'N';
  double* null = 0;

  // load M (entire matrix) into matf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> matf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = M(row, col);
      matf[index] = z.real();
      matf[index + 1] = z.imag();
    }

  zgeev_(&jobvl, &jobvr, &n, &matf[0], &n, &lambda[0], null, &n, null, &n,
         &work[0], &lwork, &rwork[0], &info);

  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in diagonalize"));
  } else if (info > 0) {
    throw(std::invalid_argument(" no convergence in diagonalize"));
  }
  eigvals.resize(n);
  for (int k = 0; k < n; k++)
    eigvals[k] = complex<double>(lambda[2 * k], lambda[2 * k + 1]);
}

// *********************************************************************

//  This takes an arbitary square complex matrix "M" and
//  calculates its determinant.  The determinant is returned
//  as a vector of complex values in "DetProd".  The product
//  of these numbers is the determinant.

void DeterminantCalculator::getDeterminantAsProduct(const CMatrix& M,
                                                    Cvector& DetProd) {
  DetProd.clear();
  int n = M.size(0);
  if (M.size(1) != uint(n)) {
    throw(std::invalid_argument("calcDeterminant requires a square matrix"));
  }
  DetProd.resize(n);
  if (n == 0) {
    return;
  }
  if (n == 1) {
    DetProd[0] = M(0, 0);
    return;
  }

  // load M (entire matrix) into matf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> matf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = M(row, col);
      matf[index] = z.real();
      matf[index + 1] = z.imag();
    }

  // do an LU decomposition:   Q = P*L*U where P is permutation
  // matrix, L is lower diagonal with 1's on diagonal, and U
  // is upper triangular
  vector<int> ipiv(n);
  int info;
  zgetrf_(&n, &n, &matf[0], &n, &ipiv[0], &info);

  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in calcDeterminant"));
  } else if (info > 0) {
    DetProd.resize(1);
    DetProd[0] = complex<double>(0.0, 0.0);
    return;
  }

  // determinant is product of diagonal elements of U
  // with (-1)^S from det(P)
  for (int k = 0; k < n; k++) {
    int index = 2 * (k + n * k);
    DetProd[k] = complex<double>(matf[index], matf[index + 1]);
  }
  // put sign from permutation matrix into first element
  int res = 1;
  for (int i = 1; i <= n; i++) {
    if (ipiv[i - 1] != i)
      res = -res;
  }
  if (res < 0)
    DetProd[0] = -DetProd[0];
}

//   returns  Omega(mu,M) where  M = complex square matrix, and
//      Omega(mu,M) = det(M) / det( sqrt(mu^2 + M M^dagger) )

complex<double> DeterminantCalculator::getOmega(double mu, const CMatrix& M) {
  int n = M.size(0);
  if (M.size(1) != uint(n))
    throw(
        std::invalid_argument("DeterminantCalculator requires square matrix"));
  if (n == 0) {
    return 0.0;
  }
  double musq = mu * mu;

  // load M (entire matrix) into nummatf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> nummatf(2 * n * n);
  for (int col = 0; col < n; col++)
    for (int row = 0; row < n; row++) {
      int index = 2 * (row + n * col);
      const complex<double>& z = M(row, col);
      nummatf[index] = z.real();
      nummatf[index + 1] = z.imag();
    }

  // load mu^2+M*Mdag (lower triangle) into denmatf fortran format
  //    (column major; row index changes fastest)
  //    complex stored as real,imag contiguous in fortran
  vector<double> denmatf(2 * n * n); // put mu^2+M*M^dag
  for (int col = 0; col < n; col++)
    for (int row = col; row < n; row++) {
      int index = 2 * (row + n * col);
      double& zre = denmatf[index];
      double& zim = denmatf[index + 1];
      for (int k = 0; k < n; k++) {
        int index1 = 2 * (row + n * k);
        double ar = nummatf[index1];
        double ai = nummatf[index1 + 1];
        int index2 = 2 * (col + n * k);
        double br = nummatf[index2];
        double bi = nummatf[index2 + 1];
        zre += ar * br + ai * bi;
        zim += ai * br - ar * bi;
      }
      if (row == col)
        zre += musq;
    }

  return get_omega(nummatf, denmatf, n);
}

//   nummatf is numerator matrix "A", denmatf is Hermitian, pos def
//   denominator matrix mu^2 + A * Adagger.  Method is to compute
//   determinant of numerator as above using LU decomposition,
//   then compute determinant of numerator using Cholesky.
//   Avoid large determinants by multiplying one diagonal of
//   numerator U and denominator L at a time.

complex<double> DeterminantCalculator::get_omega(std::vector<double>& nummatf,
                                                 std::vector<double>& denmatf,
                                                 int n) {
  int info;
  char uplo = 'L';
  zpotrf_(&uplo, &n, &denmatf[0], &n,
          &info); // Cholesky of Hermitian posdef denmatf
  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in diagonalize"));
  } else if (info > 0) {
    return 0.0;
  }

  vector<int> ipiv(n);
  zgetrf_(&n, &n, &nummatf[0], &n, &ipiv[0], &info);
  if (info < 0) {
    throw(std::invalid_argument(" bad arguments in diagonalize"));
  } else if (info > 0) {
    return 0.0;
  }

  complex<double> res(1.0, 0.0);
  for (int k = 0; k < n; k++) {
    int index = 2 * (k + n * k);
    complex<double> z(nummatf[index], nummatf[index + 1]);
    res *= z / denmatf[index];
  } // Cholesky has real pos diagonals

  for (int i = 1; i <= n; i++) {
    if (ipiv[i - 1] != i)
      res = -res;
  }
  return res;
}

// ********************************************************************
