#include "doctest.h"
#include "helper.h"
#include "fit_forms.h"
#include "param_registry.h"
#include "K_matrix_info.h"
#include "K_matrix_calc.h"
#include "xml_handler.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace TestHelper;

TEST_SUITE("ParameterNameRegistry Tests") {

    TEST_CASE("Basic Parameter Registry Operations") {
        // Get singleton instance and clear it for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Register and retrieve parameters") {
            uint hash1 = registry.registerParameter("alpha");
            uint hash2 = registry.registerParameter("beta");
            uint hash3 = registry.registerParameter("gamma_extended");
            
            CHECK(hash1 != hash2);
            CHECK(hash2 != hash3);
            CHECK(hash1 != hash3);
            
            CHECK(registry.getParameterName(hash1) == "alpha");
            CHECK(registry.getParameterName(hash2) == "beta");
            CHECK(registry.getParameterName(hash3) == "gamma_extended");
            
            CHECK(registry.hasHash(hash1));
            CHECK(registry.hasHash(hash2));
            CHECK(registry.hasHash(hash3));
            CHECK(!registry.hasHash(999999));
            
            CHECK(registry.size() == 3);
        }
        
        SUBCASE("Duplicate parameter names return same hash") {
            uint hash1 = registry.registerParameter("duplicate");
            uint hash2 = registry.registerParameter("duplicate");
            
            CHECK(hash1 == hash2);
            CHECK(registry.size() == 1);
            CHECK(registry.getParameterName(hash1) == "duplicate");
        }
        
        SUBCASE("Hash collision handling") {
            // This is difficult to test directly since we can't easily force collisions
            // but we can at least verify different names get different hashes
            std::vector<std::string> names = {"a", "b", "c", "test", "param", "x1", "y2"};
            std::vector<uint> hashes;
            
            for (const auto& name : names) {
                uint hash = registry.registerParameter(name);
                hashes.push_back(hash);
                CHECK(registry.getParameterName(hash) == name);
            }
            
            // Check all hashes are unique
            for (size_t i = 0; i < hashes.size(); ++i) {
                for (size_t j = i + 1; j < hashes.size(); ++j) {
                    CHECK(hashes[i] != hashes[j]);
                }
            }
        }
        
        SUBCASE("Non-existent hash returns empty string") {
            CHECK(registry.getParameterName(999999) == "");
        }
    }

    TEST_CASE("Expression Fit Form Integration") {
        // Clear registry for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Create Expression with parameters") {
            // Create an Expression fit form with simple parameters
            Expression expr("a_extended*x + b_extended^2");
            
            // Test parameter name detection (this works correctly)
            std::vector<std::string> paramNames = expr.getParameterNames();
            CHECK(paramNames.size() == 2);
            CHECK(paramNames[0] == "a_extended");
            CHECK(paramNames[1] == "b_extended");
            
            // Note: Parameter evaluation requires Kinitialize to be called,
            // which is private and only accessible to K-matrix calculators
        }
        
        SUBCASE("Expression with mathematical functions only") {
            Expression expr("sin(x) + cos(x)");
            
            // Test evaluation with mathematical functions (works without parameter initialization)
            std::vector<double> params; // No custom parameters
            double x = 0.0; // sin(0)=0, cos(0)=1
            double result = expr.evaluate(params, x);
            
            // Expected: sin(0) + cos(0) = 0 + 1 = 1.0
            CHECK(std::abs(result - 1.0) < 1e-10);
        }
    }

    TEST_CASE("KFitParamInfo with Expression Parameters") {
        // Clear registry for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        // Create a KElementInfo for testing
        KElementInfo kelem(0, 0, 0, 0, 0, 0, 0); // J=0, simple s-wave case
        
        SUBCASE("Create KFitParamInfo for expression parameter") {
            KFitParamInfo paramInfo(kelem, "alpha", registry);
            
            // Verify parameter was registered
            CHECK(registry.size() == 1);
            CHECK(registry.hasHash(registry.registerParameter("alpha")));
            
            // Get MCObsName and verify it can be parsed back
            std::string mcObsName = paramInfo.getMCObsName();
            CHECK(!mcObsName.empty());

            std::string retrievedName = registry.getParameterNameFromMCObsName(mcObsName);

            CHECK(retrievedName == "alpha");
        }
        
        SUBCASE("Multiple expression parameters") {
            std::vector<std::string> paramNames = {"param1", "param2", "param3"};
            std::vector<KFitParamInfo> paramInfos;
            std::vector<std::string> mcObsNames;
            
            // Create multiple expression parameters (they auto-register)
            for (const auto& name : paramNames) {
                paramInfos.emplace_back(kelem, name, registry);
                mcObsNames.push_back(paramInfos.back().getMCObsName());
            }
            
            CHECK(registry.size() == 3);
            
            // Debug: Check what hashes are actually stored in the registry
            for (uint i = 0; i < paramNames.size(); ++i) {
                // Parse the hash from the MCObsName
                uint parsedHash = registry.getHashFromMCObsName(mcObsNames[i]);
                
                // Get the hash that should be for this parameter name
                uint registeredHash = registry.registerParameter(paramNames[i]); // This should return existing hash
                
                INFO("Parameter: " << paramNames[i] 
                     << ", MCObsName: " << mcObsNames[i]
                     << ", ParsedHash: " << parsedHash 
                     << ", RegisteredHash: " << registeredHash
                     << ", HasParsedHash: " << registry.hasHash(parsedHash)
                     << ", HasRegisteredHash: " << registry.hasHash(registeredHash));
                
                std::string retrievedName = registry.getParameterNameFromMCObsName(mcObsNames[i]);
                CHECK(retrievedName == paramNames[i]);
            }
        }
    }

    TEST_CASE("Non-Expression Fit Forms Return Empty String") {
        // Clear registry for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        // Create a KElementInfo for testing
        KElementInfo kelem(0, 0, 0, 0, 0, 0, 0); // J=0, simple s-wave case
        
        SUBCASE("Polynomial parameter MCObsName") {
            KFitParamInfo polyParam(kelem, 2); // polynomial with power 2
            std::string mcObsName = polyParam.getMCObsName();
            
            // Should return empty string since it's not an expression parameter
            std::string retrievedName = registry.getParameterNameFromMCObsName(mcObsName);
            CHECK(retrievedName.empty());
        }
        
        SUBCASE("Pole energy parameter MCObsName") {
            KFitParamInfo poleParam(0, 0); // pole index 0, J=0
            std::string mcObsName = poleParam.getMCObsName();
            
            // Should return empty string since it's not an expression parameter
            std::string retrievedName = registry.getParameterNameFromMCObsName(mcObsName);
            CHECK(retrievedName.empty());
        }
        
        SUBCASE("Pole coupling parameter MCObsName") {
            KIndex kindex(0, 0, 0); // L=0, S=0, channel=0
            KFitParamInfo couplingParam(kindex, 0, 0); // pole index 0, J=0
            std::string mcObsName = couplingParam.getMCObsName();
            
            // Should return empty string since it's not an expression parameter
            std::string retrievedName = registry.getParameterNameFromMCObsName(mcObsName);
            CHECK(retrievedName.empty());
        }
    }

    TEST_CASE("MCObsName Parsing Edge Cases") {
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Invalid MCObsName formats") {
            CHECK(registry.getParameterNameFromMCObsName("").empty());
            CHECK(registry.getParameterNameFromMCObsName("NotAKStrExpr").empty());
            CHECK(registry.getParameterNameFromMCObsName("KStrExpr").empty());
            CHECK(registry.getParameterNameFromMCObsName("KStrExpr()").empty());
            CHECK(registry.getParameterNameFromMCObsName("KPoly(1,0)(0,0,0)(0,0,0)").empty());
        }
        
        SUBCASE("Valid KStrExpr but hash not in registry") {
            // This should parse correctly but return empty since hash not registered
            CHECK(registry.getParameterNameFromMCObsName("KStrExpr(999999,0)(0,0,0)(0,0,0)").empty());
        }
        
        SUBCASE("Valid KStrExpr with registered hash") {
            // Register a parameter and create a proper KStrExpr name
            uint hash = registry.registerParameter("testParam");
            
            // Create a properly formatted KStrExpr name
            std::string kstrExprName = "KStrExpr(" + std::to_string(hash) + ",0)(0,0,0)(0,0,0)";
            
            std::string retrievedName = registry.getParameterNameFromMCObsName(kstrExprName);
            CHECK(retrievedName == "testParam");
        }
    }

    TEST_CASE("Expression XML Construction") {
        // Clear registry for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Create Expression from XML") {
            // Create XML for an expression
            XMLHandler xmlin;
            xmlin.set_root("Expression");
            xmlin.put_child("String", "x + 1");
            
            Expression expr(xmlin);
            
            // Test evaluation to verify it works
            std::vector<double> params; // No custom parameters
            double x = 2.0;
            double result = expr.evaluate(params, x);
            
            // Expected: x + 1 = 2 + 1 = 3.0
            CHECK(std::abs(result - 3.0) < 1e-10);
        }
        
        SUBCASE("Expression XML output") {
            Expression expr("param1*x + param2");
            
            XMLHandler xmlout;
            expr.output(xmlout);
            
            std::string xmlString = xmlout.str();
            CHECK(xmlString.find("<Expression>") != std::string::npos);
            CHECK(xmlString.find("<String>param1*x + param2</String>") != std::string::npos);
            CHECK(xmlString.find("</Expression>") != std::string::npos);
        }
    }

    TEST_CASE("Registry Persistence") {
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Registry maintains state between accesses") {
            // Register parameters in one scope
            uint hash1 = registry.registerParameter("persistent1");
            uint hash2 = registry.registerParameter("persistent2");
            
            // Access registry again - should maintain the same state
            ParameterNameRegistry& registry2 = ParameterNameRegistry::getInstance();
            
            CHECK(&registry == &registry2); // Should be the same instance
            CHECK(registry2.size() == 2);
            CHECK(registry2.getParameterName(hash1) == "persistent1");
            CHECK(registry2.getParameterName(hash2) == "persistent2");
        }
    }

    TEST_CASE("Expression Parameter Handling Validation") {
        SUBCASE("Parameter name detection works") {
            Expression expr("a + b");
            
            // Test that parameter names are correctly identified
            std::vector<std::string> paramNames = expr.getParameterNames();
            CHECK(paramNames.size() == 2);
            CHECK(std::find(paramNames.begin(), paramNames.end(), "a") != paramNames.end());
            CHECK(std::find(paramNames.begin(), paramNames.end(), "b") != paramNames.end());
        }
        
        SUBCASE("Mathematical functions work in expressions") {
            Expression expr("sin(x) + cos(x)");
            
            // Test evaluation at x=0 where sin(0)=0, cos(0)=1
            std::vector<double> params; // No custom parameters
            double result = expr.evaluate(params, 0.0);
            
            // Expected: sin(0) + cos(0) = 0 + 1 = 1.0
            CHECK(std::abs(result - 1.0) < 1e-10);
        }
        
        SUBCASE("Complex expressions with x variable") {
            Expression expr("x^3 - 2*x^2 + x - 1");
            
            // Test with x variable (this works without parameter initialization)
            std::vector<double> params; // No custom parameters
            double x = 2.0; // 8 - 8 + 2 - 1 = 1
            double result = expr.evaluate(params, x);
            
            CHECK(std::abs(result - 1.0) < 1e-10);
        }
    }

    TEST_CASE("muParser Integration Verification") {
        SUBCASE("Simple arithmetic expression") {
            Expression expr("2 + 3");
            std::vector<double> params; // No parameters
            double result = expr.evaluate(params, 0.0);
            
            // Should get 2 + 3 = 5.0
            CHECK(std::abs(result - 5.0) < 1e-10);
        }
        
        SUBCASE("Expression with x variable") {
            Expression expr("x * 2");
            std::vector<double> params; // No parameters
            double result = expr.evaluate(params, 3.0); // x = 3.0
            
            // Should get 3.0 * 2 = 6.0
            CHECK(std::abs(result - 6.0) < 1e-10);
        }
        
        SUBCASE("Mathematical functions work") {
            Expression expr("sin(x) + cos(x)");
            std::vector<double> params; // No parameters
            double result = expr.evaluate(params, 0.0); // x = 0.0
            
            // Should get sin(0) + cos(0) = 0 + 1 = 1.0
            CHECK(std::abs(result - 1.0) < 1e-10);
        }
        
        SUBCASE("Polynomial") {
            Expression expr("x^2 + 2*x + 1");
            std::vector<double> params; // No parameters
            double result = expr.evaluate(params, 3.0); // x = 3.0
            
            // Should get 3^2 + 2*3 + 1 = 9 + 6 + 1 = 16.0
            CHECK(std::abs(result - 16.0) < 1e-10);
        }
    }

    TEST_CASE("Expression with Proper Parameter Initialization") {
        // Clear registry for clean testing
        ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
        registry.clear();
        
        SUBCASE("Expression parameter names detection") {
            // Test that Expression can correctly identify parameter names
            Expression expr1("a*x + b");
            std::vector<std::string> params1 = expr1.getParameterNames();
            CHECK(params1.size() == 2);
            CHECK(std::find(params1.begin(), params1.end(), "a") != params1.end());
            CHECK(std::find(params1.begin(), params1.end(), "b") != params1.end());
            
            Expression expr2("alpha + beta*x^2 + gamma");
            std::vector<std::string> params2 = expr2.getParameterNames();
            CHECK(params2.size() == 3);
            CHECK(std::find(params2.begin(), params2.end(), "alpha") != params2.end());
            CHECK(std::find(params2.begin(), params2.end(), "beta") != params2.end());
            CHECK(std::find(params2.begin(), params2.end(), "gamma") != params2.end());
        }
        
        SUBCASE("Expression with mathematical functions (no custom parameters)") {
            // Test expressions that work without parameter initialization
            Expression expr1("sin(x) + cos(x)");
            std::vector<double> params;
            double result1 = expr1.evaluate(params, 0.0); // sin(0) + cos(0) = 1
            CHECK(std::abs(result1 - 1.0) < 1e-10);
            
            Expression expr2("x^2 + 2*x + 1");
            double result2 = expr2.evaluate(params, 3.0); // 9 + 6 + 1 = 16
            CHECK(std::abs(result2 - 16.0) < 1e-10);
            
            Expression expr3("exp(0) + log(1)");
            double result3 = expr3.evaluate(params, 0.0); // 1 + 0 = 1
            CHECK(std::abs(result3 - 1.0) < 1e-10);
        }
        
        SUBCASE("Parameter registry interaction") {
            // Test that Expression parameters get registered correctly
            ParameterNameRegistry& registry = ParameterNameRegistry::getInstance();
            registry.clear();
            
            Expression expr("param1*x + param2");
            std::vector<std::string> paramNames = expr.getParameterNames();
            
            // Register the parameters in the registry as they would be during K-matrix setup
            for (const std::string& name : paramNames) {
                registry.registerParameter(name);
            }
            
            CHECK(registry.size() == 2);
            CHECK(registry.getParameterName(registry.registerParameter("param1")) == "param1");
            CHECK(registry.getParameterName(registry.registerParameter("param2")) == "param2");
        }
    }
} 