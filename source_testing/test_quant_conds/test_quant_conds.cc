#include "K_matrix_calc.h"
#include "box_quant.h"
#include "cnpy.h"
#include "doctest.h"
#include "helper.h"
#include "xml_handler.h"

#include <filesystem>
#include <fstream>
#include <memory>

using namespace std;
namespace fs = std::filesystem;

std::vector<double> param_vector;

// Set file_dir relative to project root
std::string file_dir =
    std::string(PROJECT_ROOT_DIR) + "/source_testing/test_quant_conds/";

const double Elab_over_mref = 4.5;

int num_of_params;

// create BoxQuantization
struct QuantHandles {
  std::unique_ptr<KtildeMatrixCalculator> kmat;
  std::unique_ptr<KtildeInverseCalculator> kinv;
  std::unique_ptr<BoxQuantization> box;
  std::unique_ptr<BoxQuantization> box_inv;
};

QuantHandles createBoxQuantization(const std::string& filename) {
  XMLHandler xml;
  xml.set_from_file(file_dir + filename);

  auto kmat = std::make_unique<KtildeMatrixCalculator>(xml);
  auto kinv = std::make_unique<KtildeInverseCalculator>(xml);

  kmat->setParameterValues(param_vector);
  kinv->setParameterValues(param_vector);

  auto box_quant = std::make_unique<BoxQuantization>(xml, kmat.get());
  auto box_quant_inv = std::make_unique<BoxQuantization>(xml, kinv.get());

  return {std::move(kmat), std::move(kinv), std::move(box_quant),
          std::move(box_quant_inv)};
}

std::vector<std::complex<double>>
hermitianMatrixToCVector(const ComplexHermitianMatrix& M) {
  const uint N = M.size();
  std::vector<std::complex<double>> vec(N * N); // N×N slots
  for (uint i = 0; i < N; ++i)
    for (uint j = 0; j < N; ++j)
      vec[i * N + j] = M.get(i, j);
  return vec;
}

std::vector<double> symmetricMatrixToVector(const RealSymmetricMatrix& M) {
  const uint N = M.size();
  std::vector<double> vec(N * N);
  for (uint i = 0; i < N; ++i)
    for (uint j = 0; j < N; ++j)
      vec[i * N + j] = M.get(i, j);
  return vec;
}

const std::string filename = "test_quant_cond.xml";

const std::string params_filename = "k_params.txt";

TEST_SUITE("Quantization Conditions Tests") {
  TEST_CASE("Load in params") {
    std::ifstream file{file_dir + params_filename};
    while (!file.eof()) {
      std::string line;
      std::getline(file, line);
      param_vector.push_back(stod(line));
      num_of_params++;
    }
  }
  TEST_CASE("BoxQuantization can be constructed") {
    auto quant_handles = createBoxQuantization(filename);
    auto box_quant = std::move(quant_handles.box);
    auto box_quant_inv = std::move(quant_handles.box_inv);

    REQUIRE(box_quant); // unique_ptr not null
    REQUIRE(box_quant_inv);
  }
  TEST_CASE("Masses can be set") {
    auto quant_handles = createBoxQuantization(filename);
    auto box_quant = std::move(quant_handles.box);
    auto box_quant_inv = std::move(quant_handles.box_inv);
    box_quant->setMassesOverRef(0, 1.5, 3.0);
    box_quant_inv->setMassesOverRef(0, 1.5, 3.0);
    CHECK(box_quant->getMass1OverRef(0) == 1.5);
    CHECK(box_quant->getMass2OverRef(0) == 3.0);
    CHECK(box_quant_inv->getMass1OverRef(0) == 1.5);
    CHECK(box_quant_inv->getMass2OverRef(0) == 3.0);
  }
  TEST_CASE("Correct number of Ktilde parameters") {
    auto quant_handles = createBoxQuantization(filename);
    auto box_quant = std::move(quant_handles.box);
    auto box_quant_inv = std::move(quant_handles.box_inv);
    CHECK(box_quant->getNumberOfKtildeParameters() == num_of_params);
  }
  TEST_CASE("Ktilde parameters can be set") {
    auto quant_handles = createBoxQuantization(filename);
    auto box_quant = std::move(quant_handles.box);
    auto box_quant_inv = std::move(quant_handles.box_inv);
    const std::vector<KFitParamInfo>& kfitparams =
        box_quant->getFitParameterInfos();
    const std::vector<KFitParamInfo>& kinvfitparams =
        box_quant_inv->getFitParameterInfos();
    for (int i = 0; i < num_of_params; ++i) {
      CHECK(box_quant->getParameterValue(kfitparams[i]) == param_vector[i]);
      CHECK(box_quant_inv->getParameterValue(kinvfitparams[i]) ==
            param_vector[i]);
    }
  }

  TEST_CASE("Quantization conditions are correct") {
    auto quant_handles = createBoxQuantization(filename);
    auto box_quant = std::move(quant_handles.box);
    auto box_quant_inv = std::move(quant_handles.box_inv);
    const int basis_size = box_quant->getBasisSize();
    const int basis_size_inv = box_quant_inv->getBasisSize();
    CMatrix CTBoxMatrix(basis_size, basis_size);
    CMatrix Stilde(basis_size, basis_size);
    CMatrix Stilde_from_inv(basis_size_inv, basis_size_inv);

    box_quant->getCTBoxMatrixFromElab(Elab_over_mref, CTBoxMatrix);
    box_quant->getStildeFromElab(Elab_over_mref, Stilde);
    box_quant_inv->getStildeFromElab(Elab_over_mref, Stilde_from_inv);

    ComplexHermitianMatrix B(basis_size, basis_size);
    ComplexHermitianMatrix B_from_inv(basis_size_inv, basis_size_inv);
    RealSymmetricMatrix Ktilde(basis_size, basis_size);
    RealSymmetricMatrix Ktilde_inv(basis_size_inv, basis_size_inv);
    box_quant->getBoxMatrixFromElab(Elab_over_mref, B);
    box_quant->getKtildeFromElab(Elab_over_mref, Ktilde);
    box_quant_inv->getKtildeinvFromElab(Elab_over_mref, Ktilde_inv);
    box_quant_inv->getBoxMatrixFromElab(Elab_over_mref, B_from_inv);
    // output to file
    std::vector<std::complex<double>> B_vec = hermitianMatrixToCVector(B);
    std::vector<std::complex<double>> B_from_inv_vec =
        hermitianMatrixToCVector(B_from_inv);
    std::vector<double> Ktilde_vec = symmetricMatrixToVector(Ktilde);
    std::vector<double> Ktilde_inv_vec = symmetricMatrixToVector(Ktilde_inv);

    std::vector<size_t> shape{B.size(), B.size()};
    std::vector<size_t> shape_inv{B_from_inv.size(), B_from_inv.size()};
    cnpy::npy_save(file_dir + "B.npy", B_vec.data(), shape, "w");
    cnpy::npy_save(file_dir + "B_from_inv.npy", B_from_inv_vec.data(),
                   shape_inv, "w");
    cnpy::npy_save(file_dir + "Ktilde.npy", Ktilde_vec.data(), shape, "w");
    cnpy::npy_save(file_dir + "Ktilde_inv.npy", Ktilde_inv_vec.data(),
                   shape_inv, "w");

    // run python script
    REQUIRE(std::system(
                ("python3 " + file_dir + "test_quant_conds.py").c_str()) == 0);

    // load in python the results
    auto refCB = cnpy::npy_load(file_dir + "CB.npy");
    auto refStilde = cnpy::npy_load(file_dir + "Stilde.npy");
    auto refStildeFromInv = cnpy::npy_load(file_dir + "Stilde_inv.npy");

    std::complex<double>* refCBvec =
        reinterpret_cast<std::complex<double>*>(refCB.data<double>());
    std::complex<double>* refStildevec =
        reinterpret_cast<std::complex<double>*>(refStilde.data<double>());
    std::complex<double>* refStildeinv =
        reinterpret_cast<std::complex<double>*>(
            refStildeFromInv.data<double>());

    CMatrix CB_ref_mat(B.size(), B.size());
    CMatrix Stilde_ref_mat(B.size(), B.size());
    CMatrix Stilde_inv_ref_mat(B_from_inv.size(), B_from_inv.size());
    CB_ref_mat.putFromVector(
        std::vector<complex<double>>(refCBvec, refCBvec + B.size() * B.size()));
    Stilde_ref_mat.putFromVector(std::vector<complex<double>>(
        refStildevec, refStildevec + B.size() * B.size()));
    Stilde_inv_ref_mat.putFromVector(std::vector<complex<double>>(
        refStildeinv, refStildeinv + B_from_inv.size() * B_from_inv.size()));

    // ensure unitarity
    CHECK(TestHelper::is_unitary(CB_ref_mat));
    CHECK(TestHelper::is_unitary(Stilde_ref_mat));
    CHECK(TestHelper::is_unitary(Stilde_inv_ref_mat));
    CHECK(TestHelper::is_unitary(CTBoxMatrix));
    CHECK(TestHelper::is_unitary(Stilde));
    CHECK(TestHelper::is_unitary(Stilde_from_inv));

    // compare
    CHECK(TestHelper::approx_equal_cmatrix(CTBoxMatrix, CB_ref_mat));
    CHECK(TestHelper::approx_equal_cmatrix(Stilde, Stilde_ref_mat));
    CHECK(
        TestHelper::approx_equal_cmatrix(Stilde_from_inv, Stilde_inv_ref_mat));

    auto qc_evs_ref = cnpy::npy_load(file_dir + "qc_evs.npy");
    auto qc_evs_inv_ref = cnpy::npy_load(file_dir + "qc_evs_inv.npy");
    auto qc_omegas_ref = cnpy::npy_load(file_dir + "qc_omegas.npy");

    std::complex<double>* qc_evs_ref_vec =
        reinterpret_cast<std::complex<double>*>(qc_evs_ref.data<double>());
    std::complex<double>* qc_evs_inv_ref_vec =
        reinterpret_cast<std::complex<double>*>(qc_evs_inv_ref.data<double>());
    std::complex<double>* qc_omega_ref_vec =
        reinterpret_cast<std::complex<double>*>(qc_omegas_ref.data<double>());

    complex<double>* evs_ref;

    for (int i = 0; i < 4; ++i) {
      BoxQuantization::QuantCondType qctype =
          static_cast<BoxQuantization::QuantCondType>(i);
      CMatrix Q(1, 1);
      CMatrix Q_ref(1, 1);
      std::vector<complex<double>> evs;
      complex<double> omega;
      double mu = 5.0;
      int stride_ind;

      if (i % 2 == 0) {
        evs = box_quant->getQCEigenvaluesFromElab(Elab_over_mref, qctype);
        evs_ref = qc_evs_ref_vec;
        stride_ind = i / 2;
        omega = box_quant->getOmegaFromElab(mu, Elab_over_mref, qctype);
      } else {
        evs = box_quant_inv->getQCEigenvaluesFromElab(Elab_over_mref, qctype);
        evs_ref = qc_evs_inv_ref_vec;
        stride_ind = (i - 1) / 2;
        omega = box_quant_inv->getOmegaFromElab(mu, Elab_over_mref, qctype);
      }

      // sort evs by increasing order in their real part
      std::sort(evs.begin(), evs.end(),
                [](const complex<double>& a, const complex<double>& b) {
                  return a.real() < b.real();
                });

      for (int j = 0; j < evs.size(); ++j) {
        string message =
            "ev " + to_string(j) + " failed for quant cond " + to_string(i);
        CHECK_MESSAGE(TestHelper::approx_equal_complex(
                          evs_ref[stride_ind * evs.size() + j], evs[j]),
                      message);
      }
      string message = "omega failed for quant cond " + to_string(i);
      CHECK_MESSAGE(
          TestHelper::approx_equal_complex(qc_omega_ref_vec[i], omega),
          message);
    }
  }
} // TEST_SUITE