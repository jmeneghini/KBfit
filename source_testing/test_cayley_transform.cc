#include "doctest.h"
#include "utils.h" // Assumed to declare calcCayleyTransformMatrix and matrix classes
#include <complex>
#include <vector>
#include <cmath>    // For std::abs, std::conj

// Helper for comparing complex numbers with tolerance
bool approx_equal_complex(const std::complex<double>& a, const std::complex<double>& b, double epsilon = 1e-6) {
    return std::abs(a.real() - b.real()) < epsilon && std::abs(a.imag() - b.imag()) < epsilon;
}

// Helper for comparing CMatrix objects
// This assumes CMatrix has size(0) for rows, size(1) for columns, and operator()(r,c) or get(r,c)
bool approx_equal_cmatrix(const CMatrix& A, const CMatrix& B, double epsilon = 1e-6) {
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

// Helper for checking if a matrix is unitary (M * M_dagger = I)
// Assumes M.size(0) is rows, M.size(1) is columns
// Assumes M.put(r,c,val) and M.get(r,c) or M(r,c)
bool is_unitary(const CMatrix& M, double epsilon = 1e-6) {
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


TEST_SUITE("calcCayleyTransformMatrix Tests") {

    TEST_CASE("ComplexHermitianMatrix: 2x2") {
        ComplexHermitianMatrix H_input(2); // Renamed from H_I for clarity
        H_input.put(0, 0, std::complex<double>(0.147, 0.0));
        H_input.put(0, 1, std::complex<double>(1.23, 2.85)); // Corrected: added closing parenthesis
        H_input.put(1, 1, std::complex<double>(1.84, 0.0));

        CMatrix CT_result;
        CMatrix CT_expected_Mform_true(2,2);
        CMatrix CT_expected_Mform_false(2,2);

        // MATHEMATICA
        // form -d*(d+i*H)*(-d+i*H)^(-1)

        // Expected result for Mform = true (d=1)
        CT_expected_Mform_true.put(0,0, std::complex<double>(-0.748231, -0.306778));
        CT_expected_Mform_true.put(0,1, std::complex<double>(-0.574324, 0.127239));
        CT_expected_Mform_true.put(1,0, std::complex<double>(0.486552 , 0.330613)); // Note: Cayley transform of Hermitian is Unitary, not necessarily Hermitian
        CT_expected_Mform_true.put(1,1, std::complex<double>(-0.808636 , 0.00832031));

        // Expected result for Mform = false (d=-1)
        CT_expected_Mform_false.put(0,0, std::complex<double>(0.748231, -0.306778));
        CT_expected_Mform_false.put(0,1, std::complex<double>(-0.486552 , 0.330613));
        CT_expected_Mform_false.put(1,0, std::complex<double>(0.574324, 0.127239));
        CT_expected_Mform_false.put(1,1, std::complex<double>(0.808636 , 0.00832031));


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
        CMatrix CT_expected_Mform_true(3,3);
        CMatrix CT_expected_Mform_false(3,3);

        // MATHEMATICA
        // form -d*(d+i*H)*(-d+i*H)^(-1)

        CT_expected_Mform_true.put(0,0, std::complex<double>(0.31383, -0.255946));
        CT_expected_Mform_true.put(0,1, std::complex<double>(-0.858359, 0.196956));
        CT_expected_Mform_true.put(0,2, std::complex<double>(0.143916, -0.199294));
        CT_expected_Mform_true.put(1,0, std::complex<double>(-0.792102, 0.0709451));
        CT_expected_Mform_true.put(1,1, std::complex<double>(-0.448677, -0.132488));
        CT_expected_Mform_true.put(1,2, std::complex<double>(-0.291093, 0.252868));
        CT_expected_Mform_true.put(2,0, std::complex<double>(0.373036, 0.253748));
        CT_expected_Mform_true.put(2,1, std::complex<double>(-0.0743473, 0.00601061));
        CT_expected_Mform_true.put(2,2, std::complex<double>(-0.88772, -0.0533493));

        CT_expected_Mform_false.put(0,0, std::complex<double>(-0.31383, -0.255946));
        CT_expected_Mform_false.put(0,1, std::complex<double>(0.792102, 0.0709451));
        CT_expected_Mform_false.put(0,2, std::complex<double>(-0.373036, 0.253748));
        CT_expected_Mform_false.put(1,0, std::complex<double>(0.858359, 0.196956));
        CT_expected_Mform_false.put(1,1, std::complex<double>(0.448677, -0.132488));
        CT_expected_Mform_false.put(1,2, std::complex<double>(0.0743473, 0.00601061));
        CT_expected_Mform_false.put(2,0, std::complex<double>(-0.143916, -0.199294));
        CT_expected_Mform_false.put(2,1, std::complex<double>(0.291093, 0.252868));
        CT_expected_Mform_false.put(2,2, std::complex<double>(0.88772, -0.0533493));

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