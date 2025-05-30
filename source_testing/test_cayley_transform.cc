#include "doctest.h"
#include "helper.h"
#include "utils.h" // Assumed to declare calcCayleyTransformMatrix and matrix classes
#include <cmath> // For std::abs, std::conj
#include <complex>
#include <vector>

using namespace TestHelper;

TEST_SUITE("calcCayleyTransformMatrix Tests") {

  TEST_CASE("ComplexHermitianMatrix: 2x2") {
    ComplexHermitianMatrix H_input(2);
    H_input.put(0, 0, std::complex<double>(0.147, 0.0));
    H_input.put(0, 1, std::complex<double>(1.23, 2.85));
    H_input.put(1, 1, std::complex<double>(1.84, 0.0));

    CMatrix CT_result;
    CMatrix CT_expected_Mform_true(2, 2);
    CMatrix CT_expected_Mform_false(2, 2);

    // MATHEMATICA
    // form -d*(d+i*H)*(-d+i*H)^(-1)

    // Expected result for Mform = true (d=1)
    CT_expected_Mform_true.put(0, 0,
                               std::complex<double>(-0.748231, -0.306778));
    CT_expected_Mform_true.put(0, 1, std::complex<double>(-0.574324, 0.127239));
    CT_expected_Mform_true.put(1, 0, std::complex<double>(0.486552, 0.330613));
    CT_expected_Mform_true.put(1, 1,
                               std::complex<double>(-0.808636, 0.00832031));

    // Expected result for Mform = false (d=-1)
    CT_expected_Mform_false.put(0, 0,
                                std::complex<double>(0.748231, -0.306778));
    CT_expected_Mform_false.put(0, 1,
                                std::complex<double>(-0.486552, 0.330613));
    CT_expected_Mform_false.put(1, 0, std::complex<double>(0.574324, 0.127239));
    CT_expected_Mform_false.put(1, 1,
                                std::complex<double>(0.808636, 0.00832031));

    SUBCASE("Mform = true (d=1)") {
      calcCayleyTransformMatrix(H_input, CT_result, true);
      // The Cayley transform of a Hermitian matrix should be unitary.
      CHECK(is_unitary(CT_result));
      CHECK(approx_equal_cmatrix(CT_result, CT_expected_Mform_true));
    }

    SUBCASE("Mform = false (d=-1)") {
      calcCayleyTransformMatrix(H_input, CT_result, false);
      // The Cayley transform of a Hermitian matrix should be unitary.
      CHECK(is_unitary(CT_result));
      CHECK(approx_equal_cmatrix(CT_result, CT_expected_Mform_false));
    }
  }

  TEST_CASE("ComplexHermitianMatrix: 3x3") {
    ComplexHermitianMatrix H_input(3); // Renamed from H_I for clarity
    H_input.put(0, 0, std::complex<double>(0.5, 0.0));
    H_input.put(1, 1, std::complex<double>(1.5, 0.0));
    H_input.put(2, 2, std::complex<double>(2.5, 0.0));

    H_input.put(0, 1, std::complex<double>(1.1, 2.2));
    H_input.put(0, 2, std::complex<double>(3.3, 4.4));
    H_input.put(1, 2, std::complex<double>(5.5, 6.6));

    CMatrix CT_result;
    CMatrix CT_expected_Mform_true(3, 3);
    CMatrix CT_expected_Mform_false(3, 3);

    // MATHEMATICA
    // form -d*(d+i*H)*(-d+i*H)^(-1)

    CT_expected_Mform_true.put(0, 0, std::complex<double>(0.31383, -0.255946));
    CT_expected_Mform_true.put(0, 1, std::complex<double>(-0.858359, 0.196956));
    CT_expected_Mform_true.put(0, 2, std::complex<double>(0.143916, -0.199294));
    CT_expected_Mform_true.put(1, 0,
                               std::complex<double>(-0.792102, 0.0709451));
    CT_expected_Mform_true.put(1, 1,
                               std::complex<double>(-0.448677, -0.132488));
    CT_expected_Mform_true.put(1, 2, std::complex<double>(-0.291093, 0.252868));
    CT_expected_Mform_true.put(2, 0, std::complex<double>(0.373036, 0.253748));
    CT_expected_Mform_true.put(2, 1,
                               std::complex<double>(-0.0743473, 0.00601061));
    CT_expected_Mform_true.put(2, 2,
                               std::complex<double>(-0.88772, -0.0533493));

    CT_expected_Mform_false.put(0, 0,
                                std::complex<double>(-0.31383, -0.255946));
    CT_expected_Mform_false.put(0, 1,
                                std::complex<double>(0.792102, 0.0709451));
    CT_expected_Mform_false.put(0, 2,
                                std::complex<double>(-0.373036, 0.253748));
    CT_expected_Mform_false.put(1, 0, std::complex<double>(0.858359, 0.196956));
    CT_expected_Mform_false.put(1, 1,
                                std::complex<double>(0.448677, -0.132488));
    CT_expected_Mform_false.put(1, 2,
                                std::complex<double>(0.0743473, 0.00601061));
    CT_expected_Mform_false.put(2, 0,
                                std::complex<double>(-0.143916, -0.199294));
    CT_expected_Mform_false.put(2, 1, std::complex<double>(0.291093, 0.252868));
    CT_expected_Mform_false.put(2, 2,
                                std::complex<double>(0.88772, -0.0533493));

    SUBCASE("Mform = true (d=1)") {
      calcCayleyTransformMatrix(H_input, CT_result, true);
      // The Cayley transform of a Hermitian matrix should be unitary.
      CHECK(is_unitary(CT_result));
      CHECK(approx_equal_cmatrix(CT_result, CT_expected_Mform_true));
    }

    SUBCASE("Mform = false (d=-1)") {
      calcCayleyTransformMatrix(H_input, CT_result, false);
      // The Cayley transform of a Hermitian matrix should be unitary.
      CHECK(is_unitary(CT_result));
      CHECK(approx_equal_cmatrix(CT_result, CT_expected_Mform_false));
    }
  }
}