#include "matrix.h"
#include "utils.h"
#include "xml_handler.h"
#include <random>

using namespace std;

complex<double> get_random_complex(bool imagzero) {
  if (imagzero)
    return complex<double>(double(rand() % 1048576) / 1048576.0 - 0.5, 0.0);
  else
    return complex<double>(double(rand() % 1048576) / 1048576.0 - 0.5,
                           double(rand() % 1048576) / 1048576.0 - 0.5);
}

void get_random_herm_matrix(ComplexHermitianMatrix& H, uint n) {
  H.resize(n);
  for (uint row = 0; row < n; ++row) {
    H.put(row, row, get_random_complex(true));
    for (uint col = 0; col < row; ++col)
      H.put(row, col, get_random_complex(false));
  }
}

// ***************************************************************

bool test_a_matrix_inverse(const ComplexHermitianMatrix& H) {
  double eps = 1e-10;
  uint n = H.size();
  cout << endl << "Dimension of matrix = " << n << endl;
  ComplexHermitianMatrix Hinv;
  calcMatrixInverse(H, Hinv);
  // for (uint row=0;row<n;++row)
  // for (uint col=0;col<=row;++col)
  //    cout << "H("<<row<<","<<col<<") = "<<H(row,col)<<endl;
  // for (uint row=0;row<n;++row)
  // for (uint col=0;col<=row;++col)
  //    cout << "Hinv("<<row<<","<<col<<") = "<<Hinv(row,col)<<endl;
  bool flag = true;
  complex<double> zero(0.0, 0.0);
  complex<double> one(1.0, 0.0);
  complex<double> diff;
  // check H*Hinv=1
  for (uint row = 0; row < n; ++row)
    for (uint col = 0; col < n; ++col) {
      complex<double> res(0.0, 0.0);
      for (uint k = 0; k < n; ++k)
        res += H(row, k) * Hinv(k, col);
      complex<double>& correct = (row == col) ? one : zero;
      diff = res - correct;
      if (abs(diff) > eps) {
        cout << "diff = " << diff << endl;
        flag = false;
      }
    }
  // check Hinv*H=1
  for (uint row = 0; row < n; ++row)
    for (uint col = 0; col < n; ++col) {
      complex<double> res(0.0, 0.0);
      for (uint k = 0; k < n; ++k)
        res += Hinv(row, k) * H(k, col);
      complex<double>& correct = (row == col) ? one : zero;
      diff = res - correct;
      if (abs(diff) > eps)
        flag = false;
    }
  if (flag)
    cout << "SUCCESS" << endl << endl;
  else
    cout << "FAILURE" << endl << endl;
  return flag;
}

// ***************************************************************

void testMatrixInverse(XMLHandler& xml_in) {
  if (xml_tag_count(xml_in, "TestMatrixInverse") == 0)
    return;

  cout << "TestMatrixInverse" << endl;
  cout.precision(12);

  cout << "One-dimensional test:" << endl;
  ComplexHermitianMatrix H(1);
  H.put(0, 0, complex<double>(3.2, 0.0));
  bool flag = test_a_matrix_inverse(H);

  cout << "Two-dimensional test:" << endl;
  H.resize(2);
  H.put(0, 0, complex<double>(4.1, 0.0));
  H.put(1, 1, complex<double>(-2.3, 0.0));
  H.put(0, 1, complex<double>(-1.7, 2.9));
  flag &= test_a_matrix_inverse(H);

  cout << "Three-dimensional test:" << endl;
  H.resize(3);
  H.put(0, 0, complex<double>(4.1, 0.0));
  H.put(1, 1, complex<double>(-2.3, 0.0));
  H.put(2, 2, complex<double>(5.5, 0.0));
  H.put(0, 1, complex<double>(3.7, -4.9));
  H.put(0, 2, complex<double>(-6.6, 1.9));
  H.put(1, 2, complex<double>(8.1, 3.3));
  flag &= test_a_matrix_inverse(H);

  for (uint n = 4; n < 36; n += 5) {
    get_random_herm_matrix(H, n);
    flag &= test_a_matrix_inverse(H);
  }

  if (flag)
    cout << "ALL SUCCESS" << endl;
  else
    cout << "At least one FAILURE" << endl;
}

// ***************************************************************
