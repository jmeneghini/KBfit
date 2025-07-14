#include "K_matrix_calc.h"
#include "box_quant.h"
#include "chisq_detres.h"
#include "chisq_spectrum.h"
#include "doctest.h"
#include "helper.h"
#include "xml_handler.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace TestHelper;
namespace fs = std::filesystem;

namespace {
// Set file_dir relative to project root
std::string file_dir = std::string(PROJECT_ROOT_DIR) + "/source_testing/test_quant_conds/";
}

std::vector<double> loadParametersFromFile(const std::string& filename) {
    std::vector<double> params;
    std::ifstream file(file_dir + filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open parameter file: " + filename);
    }
    
    double value;
    while (file >> value) {
        params.push_back(value);
    }
    file.close();
    return params;
}

struct CloneTestHandles {
    std::unique_ptr<KtildeMatrixCalculator> kmat_original;
    std::unique_ptr<KtildeInverseCalculator> kinv_original;
    std::unique_ptr<BoxQuantization> box_original;
    std::unique_ptr<BoxQuantization> box_inv_original;
    
    std::unique_ptr<KtildeMatrixCalculator> kmat_cloned;
    std::unique_ptr<KtildeInverseCalculator> kinv_cloned;
    std::unique_ptr<BoxQuantization> box_cloned;
    std::unique_ptr<BoxQuantization> box_inv_cloned;
};

CloneTestHandles createCloneTestHandles(const std::string& xml_filename, 
                                      const std::vector<double>& params) {
    XMLHandler xml;
    xml.set_from_file(file_dir + xml_filename);

    CloneTestHandles handles;
    
    // Create original objects
    handles.kmat_original = std::make_unique<KtildeMatrixCalculator>(xml);
    handles.kinv_original = std::make_unique<KtildeInverseCalculator>(xml);
    
    handles.kmat_original->setParameterValues(params);
    handles.kinv_original->setParameterValues(params);
    
    handles.box_original = std::make_unique<BoxQuantization>(xml, handles.kmat_original.get());
    handles.box_inv_original = std::make_unique<BoxQuantization>(xml, handles.kinv_original.get());
    
    // Create cloned objects
    handles.kmat_cloned = handles.kmat_original->clone();
    handles.kinv_cloned = handles.kinv_original->clone();
    
    handles.box_cloned = handles.box_original->clone(handles.kmat_cloned.get());
    handles.box_inv_cloned = handles.box_inv_original->clone(nullptr, handles.kinv_cloned.get());
    
    return handles;
}

TEST_SUITE("Clone Methods Tests") {

    TEST_CASE("KtildeMatrixCalculator Clone Tests with Real XML") {
        auto params = loadParametersFromFile("k_params.txt");
        
        SUBCASE("Basic clone properties") {
            XMLHandler xml;
            xml.set_from_file(file_dir + "test_quant_cond.xml");
            
            KtildeMatrixCalculator original(xml);
            original.setParameterValues(params);
            
            auto cloned = original.clone();
            
            // Test that objects are different
            CHECK(cloned.get() != &original);
            
            // Test that properties are identical
            CHECK(cloned->getNumberOfParameters() == original.getNumberOfParameters());
            CHECK(cloned->getNumberOfDecayChannels() == original.getNumberOfDecayChannels());
            CHECK(cloned->getElementInfos() == original.getElementInfos());
            CHECK(cloned->getFitParameterInfos() == original.getFitParameterInfos());
            
            // Test parameter independence
            std::vector<double> new_params = params;
            for (auto& p : new_params) p *= 1.5;
            
            cloned->setParameterValues(new_params);
            CHECK(original.getParameterValues() != cloned->getParameterValues());
            CHECK(cloned->getParameterValues() == new_params);
            
            // Test that calculation results are different with different parameters
            double Ecm = 2.5;
            double orig_result = original.calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            double clone_result = cloned->calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            
            // Should be different since we changed parameters
            CHECK(std::abs(orig_result - clone_result) > 1e-10);
        }
        
        SUBCASE("Clone with same parameters produces same results") {
            XMLHandler xml;
            xml.set_from_file(file_dir + "test_quant_cond.xml");
            
            KtildeMatrixCalculator original(xml);
            original.setParameterValues(params);
            
            auto cloned = original.clone();
            cloned->setParameterValues(params); // Set same parameters
            
            double Ecm = 2.5;
            double orig_result = original.calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            double clone_result = cloned->calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            
            // Should be identical with same parameters
            CHECK(std::abs(orig_result - clone_result) < 1e-15);
        }
    }

    TEST_CASE("KtildeInverseCalculator Clone Tests with Real XML") {
        auto params = loadParametersFromFile("k_params.txt");
        
        SUBCASE("Basic clone properties") {
            XMLHandler xml;
            xml.set_from_file(file_dir + "test_quant_cond.xml");
            
            KtildeInverseCalculator original(xml);
            original.setParameterValues(params);
            
            auto cloned = original.clone();
            
            // Test that objects are different
            CHECK(cloned.get() != &original);
            
            // Test that properties are identical
            CHECK(cloned->getNumberOfParameters() == original.getNumberOfParameters());
            CHECK(cloned->getNumberOfDecayChannels() == original.getNumberOfDecayChannels());
            CHECK(cloned->getElementInfos() == original.getElementInfos());
            
            // Test calculation results with same parameters
            cloned->setParameterValues(params);
            double Ecm = 2.5;
            double orig_result = original.calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            double clone_result = cloned->calculate(0, 0, 0, 0, 0, 0, 0, Ecm);
            
            CHECK(std::abs(orig_result - clone_result) < 1e-15);
        }
    }

    TEST_CASE("BoxQuantization Clone Tests with Real XML") {
        auto params = loadParametersFromFile("k_params.txt");
        
        SUBCASE("Clone with KtildeMatrix") {
            auto handles = createCloneTestHandles("test_quant_cond.xml", params);
            
            // Test basic properties
            CHECK(handles.box_cloned.get() != handles.box_original.get());
            CHECK(handles.box_cloned->getMomRay() == handles.box_original->getMomRay());
            CHECK(handles.box_cloned->getLittleGroupIrrep() == handles.box_original->getLittleGroupIrrep());
            CHECK(handles.box_cloned->getTotalMomentumIntegerSquared() == handles.box_original->getTotalMomentumIntegerSquared());
            CHECK(handles.box_cloned->getNumberOfDecayChannels() == handles.box_original->getNumberOfDecayChannels());
            CHECK(handles.box_cloned->getBasisSize() == handles.box_original->getBasisSize());
            
            // Test mode consistency
            CHECK(handles.box_cloned->isKtildeMode() == handles.box_original->isKtildeMode());
            CHECK(handles.box_cloned->isKtildeInverseMode() == handles.box_original->isKtildeInverseMode());
            
            // Test mass and length independence
            handles.box_original->setRefMassL(6.0);
            handles.box_original->setMassesOverRef(0, 1.0, 1.0);
            handles.box_original->setMassesOverRef(1, 1.2, 1.2);
            
            handles.box_cloned->setRefMassL(8.0);
            handles.box_cloned->setMassesOverRef(0, 1.5, 1.5);
            handles.box_cloned->setMassesOverRef(1, 1.8, 1.8);
            
            CHECK(std::abs(handles.box_original->getRefMassL() - 6.0) < 1e-10);
            CHECK(std::abs(handles.box_cloned->getRefMassL() - 8.0) < 1e-10);
            CHECK(std::abs(handles.box_original->getMass1OverRef(0) - 1.0) < 1e-10);
            CHECK(std::abs(handles.box_cloned->getMass1OverRef(0) - 1.5) < 1e-10);
        }
        
        SUBCASE("Clone with KtildeInverse") {
            auto handles = createCloneTestHandles("test_quant_cond.xml", params);
            
            // Test that the inverse clone works properly
            CHECK(handles.box_inv_cloned.get() != handles.box_inv_original.get());
            CHECK(handles.box_inv_cloned->isKtildeInverseMode() == handles.box_inv_original->isKtildeInverseMode());
            CHECK(handles.box_inv_cloned->getBasisSize() == handles.box_inv_original->getBasisSize());
        }
        
        SUBCASE("Matrix calculations with cloned objects") {
            auto handles = createCloneTestHandles("test_quant_cond.xml", params);
            
            // Set same parameters for fair comparison
            handles.box_original->setRefMassL(6.0);
            handles.box_original->setMassesOverRef(0, 1.0, 1.0);
            handles.box_original->setMassesOverRef(1, 1.2, 1.2);
            
            handles.box_cloned->setRefMassL(6.0);
            handles.box_cloned->setMassesOverRef(0, 1.0, 1.0);
            handles.box_cloned->setMassesOverRef(1, 1.2, 1.2);
            
            double Elab = 4.5;
            
            // Test box matrix calculation
            ComplexHermitianMatrix B_orig, B_cloned;
            handles.box_original->getBoxMatrixFromElab(Elab, B_orig);
            handles.box_cloned->getBoxMatrixFromElab(Elab, B_cloned);
            
            CHECK(B_orig.size() == B_cloned.size());
            for (uint i = 0; i < B_orig.size(); ++i) {
                for (uint j = 0; j <= i; ++j) {
                    CHECK(std::abs(B_orig(i, j) - B_cloned(i, j)) < 1e-12);
                }
            }
            
            // Test K-matrix calculation
            RealSymmetricMatrix K_orig, K_cloned;
            handles.box_original->getKtildeFromElab(Elab, K_orig);
            handles.box_cloned->getKtildeFromElab(Elab, K_cloned);
            
            CHECK(K_orig.size() == K_cloned.size());
            for (uint i = 0; i < K_orig.size(); ++i) {
                for (uint j = 0; j <= i; ++j) {
                    CHECK(std::abs(K_orig(i, j) - K_cloned(i, j)) < 1e-12);
                }
            }
        }
    }

    TEST_CASE("Integration Tests with QuantCondType") {
        auto params = loadParametersFromFile("k_params.txt");
        auto handles = createCloneTestHandles("test_quant_cond.xml", params);
        
        // Set same physical parameters
        double Elab = 4.5;
        double mu = 5.0;
        
        handles.box_original->setRefMassL(6.0);
        handles.box_original->setMassesOverRef(0, 1.0, 1.0);
        handles.box_original->setMassesOverRef(1, 1.2, 1.2);
        
        handles.box_cloned->setRefMassL(6.0);
        handles.box_cloned->setMassesOverRef(0, 1.0, 1.0);
        handles.box_cloned->setMassesOverRef(1, 1.2, 1.2);
        
        handles.box_inv_original->setRefMassL(6.0);
        handles.box_inv_original->setMassesOverRef(0, 1.0, 1.0);
        handles.box_inv_original->setMassesOverRef(1, 1.2, 1.2);
        
        handles.box_inv_cloned->setRefMassL(6.0);
        handles.box_inv_cloned->setMassesOverRef(0, 1.0, 1.0);
        handles.box_inv_cloned->setMassesOverRef(1, 1.2, 1.2);
        
        SUBCASE("Omega calculations consistency") {
            // Test KtildeB quantization condition
            auto omega_orig = handles.box_original->getOmegaFromElab(mu, Elab, BoxQuantization::KtildeB);
            auto omega_cloned = handles.box_cloned->getOmegaFromElab(mu, Elab, BoxQuantization::KtildeB);
            
            CHECK(std::abs(omega_orig - omega_cloned) < 1e-12);
            
            // Test KtildeinvB quantization condition
            auto omega_inv_orig = handles.box_inv_original->getOmegaFromElab(mu, Elab, BoxQuantization::KtildeinvB);
            auto omega_inv_cloned = handles.box_inv_cloned->getOmegaFromElab(mu, Elab, BoxQuantization::KtildeinvB);
            
            CHECK(std::abs(omega_inv_orig - omega_inv_cloned) < 1e-12);
        }
        
        SUBCASE("Eigenvalue calculations consistency") {
            auto evs_orig = handles.box_original->getQCEigenvaluesFromElab(Elab, BoxQuantization::KtildeB);
            auto evs_cloned = handles.box_cloned->getQCEigenvaluesFromElab(Elab, BoxQuantization::KtildeB);
            
            CHECK(evs_orig.size() == evs_cloned.size());
            
            // Sort both vectors for comparison
            std::sort(evs_orig.begin(), evs_orig.end(),
                      [](const complex<double>& a, const complex<double>& b) {
                          return a.real() < b.real();
                      });
            std::sort(evs_cloned.begin(), evs_cloned.end(),
                      [](const complex<double>& a, const complex<double>& b) {
                          return a.real() < b.real();
                      });
            
            for (size_t i = 0; i < evs_orig.size(); ++i) {
                CHECK(std::abs(evs_orig[i] - evs_cloned[i]) < 1e-12);
            }
        }
    }

    TEST_CASE("Memory Safety and Independence") {
        auto params = loadParametersFromFile("k_params.txt");
        
        SUBCASE("Verify deep copy behavior") {
            std::unique_ptr<KtildeMatrixCalculator> cloned;
            
            {
                // Create original in limited scope
                XMLHandler xml;
                xml.set_from_file(file_dir + "test_quant_cond.xml");
                KtildeMatrixCalculator original(xml);
                original.setParameterValues(params);
                
                cloned = original.clone();
                
                // Verify clone works while original exists
                CHECK(cloned->getNumberOfParameters() == original.getNumberOfParameters());
            } // original goes out of scope here
            
            // Clone should still work after original is destroyed
            CHECK(cloned->getNumberOfParameters() > 0);
            CHECK(cloned->getNumberOfDecayChannels() == 2);  // From test XML
            
            // Should be able to use cloned object
            double result = cloned->calculate(0, 0, 0, 0, 0, 0, 0, 2.5);
            CHECK(std::isfinite(result));  // Should produce valid result
        }
        
        SUBCASE("Parameter independence after cloning") {
            XMLHandler xml;
            xml.set_from_file(file_dir + "test_quant_cond.xml");
            
            KtildeMatrixCalculator original(xml);
            original.setParameterValues(params);
            
            auto cloned = original.clone();
            
            // Modify original parameters
            std::vector<double> new_params = params;
            for (auto& p : new_params) p *= 2.0;
            original.setParameterValues(new_params);
            
            // Clone should still have original parameters (until explicitly changed)
            cloned->setParameterValues(params);
            
            // Verify they produce different results
            double result_orig = original.calculate(0, 0, 0, 0, 0, 0, 0, 2.5);
            double result_clone = cloned->calculate(0, 0, 0, 0, 0, 0, 0, 2.5);
            
            CHECK(std::abs(result_orig - result_clone) > 1e-10);
        }
    }
} 