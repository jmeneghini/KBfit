#include "fit_forms.h"
#include <cmath>
#include <cctype>
#include <algorithm>
#include <regex>
#include <set>
#include <limits>

using namespace std;

// ****************************************************************

std::string FitForm::output(int indent) const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.output(indent);
}

std::string FitForm::str() const {
  XMLHandler xmlout;
  output(xmlout);
  return xmlout.str();
}

// ****************************************************************

Polynomial::Polynomial(uint degree) {
  for (uint k = 0; k <= degree; k++)
    m_powers.insert(k);
}

Polynomial::Polynomial(const std::set<uint>& powers) : m_powers(powers) {
  if (powers.empty())
    throw(std::invalid_argument(
        "One power of nonzero coefficient required in Polynomial"));
}

Polynomial::Polynomial(XMLHandler& xmlin) { initialize(xmlin); }

void Polynomial::initialize(XMLHandler& xmlin) {
  m_powers.clear();
  m_coef_indices.clear();
  try {
    XMLHandler xmlf(xmlin, "Polynomial");
    uint count1 = xmlf.count_among_children("Degree");
    uint count2 = xmlf.count_among_children("Powers");
    if ((count1 + count2) != 1)
      throw(std::invalid_argument(
          "Polynomial must have either <Degree> or <Powers> tag"));
    if (count1 == 1) {
      uint deg;
      xmlreadchild(xmlf, "Degree", deg);
      for (uint k = 0; k <= deg; k++)
        m_powers.insert(k);
    } else {
      vector<uint> ipowers;
      xmlreadchild(xmlf, "Powers", ipowers);
      if (ipowers.empty())
        throw(std::invalid_argument(
            "One power of nonzero coefficient required in Polynomial"));
      for (uint k = 0; k < ipowers.size(); k++)
        m_powers.insert(ipowers[k]);
    }
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not initialize Polynomial: ") +
                                xp.what()));
  }
}

void Polynomial::output(XMLHandler& xmlout) const {
  xmlout.set_root("Polynomial");
  uint deg = 0;
  for (set<uint>::iterator it = m_powers.begin(); it != m_powers.end(); it++)
    if ((*it) > deg)
      deg = (*it);
  if (m_powers.size() == (deg + 1)) {
    xmlout.put_child("Degree", make_string(deg));
    return;
  }
  vector<uint> ipow(m_powers.begin(), m_powers.end());
  xmlout.put_child("Powers", make_string(ipow));
}

//  sets up m_coef_indices and modifies the paramindices map

void Polynomial::Kinitialize(const KElementInfo& kinfo,
                             std::map<KFitParamInfo, uint>& paramindices) {
  uint deg = 0;
  for (set<uint>::iterator st = m_powers.begin(); st != m_powers.end(); st++)
    if ((*st) > deg)
      deg = (*st);
  m_coef_indices.resize(deg + 1);
  for (uint k = 0; k <= deg; k++) {
    if (m_powers.find(k) == m_powers.end())
      m_coef_indices[k] = -1;
    else {
      uint nextindex = paramindices.size();
      m_coef_indices[k] = nextindex;
      KFitParamInfo kpinfo(kinfo, k);
      paramindices.insert(make_pair(kpinfo, nextindex));
    }
  }
}

// ***************************************************************************************

SumOfPoles::SumOfPoles(uint numpoles) {
  if (numpoles < 1)
    throw(std::invalid_argument("Must have at least one pole in SumOfPoles"));
  for (uint k = 0; k < numpoles; k++)
    m_poleindices.insert(k);
}

SumOfPoles::SumOfPoles(const std::set<uint>& poleindices)
    : m_poleindices(poleindices) {
  if (poleindices.empty())
    throw(std::invalid_argument("Must have at least one pole in SumOfPoles"));
}

SumOfPoles::SumOfPoles(XMLHandler& xmlin) { initialize(xmlin); }

void SumOfPoles::initialize(XMLHandler& xmlin) {
  m_poleindices.clear();
  m_rowg.clear();
  m_colg.clear();
  m_poles.clear();
  try {
    XMLHandler xmlf(xmlin, "SumOfPoles");
    uint count1 = xmlf.count_among_children("NumberOfPoles");
    uint count2 = xmlf.count_among_children("PoleIndices");
    if ((count1 + count2) != 1)
      throw(std::invalid_argument(
          "SumOfPoles must have either <NumberOfPoles> or <PoleIndices> tag"));
    if (count1 == 1) {
      uint numpoles;
      xmlreadchild(xmlin, "NumberOfPoles", numpoles);
      if (numpoles < 1)
        throw(
            std::invalid_argument("Must have at least one pole in SumOfPoles"));
      for (uint k = 0; k < numpoles; k++)
        m_poleindices.insert(k);
    } else {
      vector<uint> ipoles;
      xmlreadchild(xmlf, "PoleIndices", ipoles);
      if (ipoles.empty())
        throw(
            std::invalid_argument("Must have at least one pole in SumOfPoles"));
      for (uint k = 0; k < ipoles.size(); k++)
        m_poleindices.insert(ipoles[k]);
    }
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(string("Could not initialize SumOfPoles: ") +
                                xp.what()));
  }
}

void SumOfPoles::output(XMLHandler& xmlout) const {
  xmlout.set_root("SumOfPoles");
  uint maxpoles = 0;
  for (set<uint>::iterator it = m_poleindices.begin();
       it != m_poleindices.end(); it++)
    if ((*it) > maxpoles)
      maxpoles = (*it);
  if (m_poleindices.size() == (maxpoles + 1)) {
    xmlout.put_child("NumberOfPoles", make_string(m_poleindices.size()));
    return;
  }
  vector<uint> ipow(m_poleindices.begin(), m_poleindices.end());
  xmlout.put_child("PoleIndices", make_string(ipow));
}

//  sets up m_rowg, m_colg, m_poles, and modifies the paramindices map

void SumOfPoles::Kinitialize(const KElementInfo& kinfo,
                             std::map<KFitParamInfo, uint>& paramindices) {
  uint npoles = m_poleindices.size();
  m_rowg.resize(npoles);
  m_colg.resize(npoles);
  m_poles.resize(npoles);
  std::map<KFitParamInfo, uint>::iterator kt;
  vector<uint> poleindices(m_poleindices.begin(), m_poleindices.end());
  for (uint k = 0; k < npoles; k++) {
    uint Jtimestwo = kinfo.getJtimestwo();
    KFitParamInfo kpenergy(poleindices[k], Jtimestwo);
    KFitParamInfo krowg(kinfo.getRow(), poleindices[k], Jtimestwo);
    KFitParamInfo kcolg(kinfo.getColumn(), poleindices[k], Jtimestwo);
    vector<KFitParamInfo*> kptrs(3);
    kptrs[0] = &kpenergy;
    kptrs[1] = &krowg;
    kptrs[2] = &kcolg;
    vector<uint*> iptrs(3);
    iptrs[0] = &(m_poles[k]);
    iptrs[1] = &(m_rowg[k]);
    iptrs[2] = &(m_colg[k]);
    for (uint kk = 0; kk < 3; kk++) {
      kt = paramindices.find(*(kptrs[kk]));
      if (kt != paramindices.end()) {
        (*(iptrs[kk])) = kt->second;
      } else {
        (*(iptrs[kk])) = paramindices.size();
        paramindices.insert(make_pair(*(kptrs[kk]), (*(iptrs[kk]))));
      }
    }
  }
}

// ***************************************************************************************

SumOfPolesPlusPolynomial::SumOfPolesPlusPolynomial(XMLHandler& xmlin) {
  try {
    XMLHandler xmlf(xmlin, "SumOfPolesPlusPolynomial");
    m_sumpoles.initialize(xmlf);
    m_poly.initialize(xmlf);
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(
        string("Could not initialize SumOfPolesPlusPolynomial: ") + xp.what()));
  }
}

void SumOfPolesPlusPolynomial::output(XMLHandler& xmlout) const {
  xmlout.set_root("SumOfPolesPlusPolynomial");
  XMLHandler xmls;
  m_sumpoles.output(xmls);
  xmlout.put_child(xmls);
  m_poly.output(xmls);
  xmlout.put_child(xmls);
}

//  sets up m_rowg, m_colg, m_poles, and modifies the paramindices map

void SumOfPolesPlusPolynomial::Kinitialize(
    const KElementInfo& kinfo, std::map<KFitParamInfo, uint>& paramindices) {
  m_sumpoles.Kinitialize(kinfo, paramindices);
  m_poly.Kinitialize(kinfo, paramindices);
}

// ***************************************************************************************

Expression::Expression() 
    : m_expression("0"), m_x_value(0.0), m_parser_initialized(false) {}

Expression::Expression(const std::string& expression)
    : m_expression(expression), m_x_value(0.0), m_parser_initialized(false) {
  parseExpression();
}

Expression::Expression(XMLHandler& xmlin) : m_parser_initialized(false) {
  try {
    XMLHandler xmlf(xmlin, "Expression");
    
    // Read expression
    xmlread(xmlf, "String", m_expression, "Expression");
    
    m_x_value = 0.0;
    
    // Parse the expression to extract parameter names
    parseExpression();
  } catch (const std::exception& xp) {
    throw(std::invalid_argument(
        std::string("Could not initialize Expression: ") + xp.what()));
  }
}

Expression::Expression(const Expression& in)
    : m_expression(in.m_expression), 
      m_param_names(in.m_param_names), m_param_indices(in.m_param_indices),
      m_x_value(0.0), m_parser_initialized(false) {}

Expression& Expression::operator=(const Expression& in) {
  m_expression = in.m_expression;
  m_param_names = in.m_param_names;
  m_param_indices = in.m_param_indices;
  m_x_value = 0.0;
  m_parser_initialized = false;
  return *this;
}

void Expression::parseExpression() {
  // Extract variable names from the expression string using regex
  // This approach finds all potential variable names before setting up muParser
  try {
    m_param_names.clear();
    
    // Use regex to find all potential variable names (alphanumeric + underscore, not starting with digit)
    std::regex var_regex(R"([a-zA-Z_][a-zA-Z0-9_]*)");
    std::sregex_iterator iter(m_expression.begin(), m_expression.end(), var_regex);
    std::sregex_iterator end;
    
    std::set<std::string> unique_vars;
    for (; iter != end; ++iter) {
      std::string var = iter->str();
      // Skip common mathematical functions and constants
      if (var != "x" && var != "sin" && var != "cos" && var != "tan" && 
          var != "exp" && var != "log" && var != "sqrt" && var != "abs" &&
          var != "pi" && var != "e" && var != "ln" && var != "log10" &&
          var != "sinh" && var != "cosh" && var != "tanh" && var != "asin" &&
          var != "acos" && var != "atan" && var != "atan2" && var != "pow" &&
          var != "min" && var != "max" && var != "sum" && var != "avg") {
        unique_vars.insert(var);
      }
    }
    
    // Convert set to vector
    for (const auto& var : unique_vars) {
      m_param_names.push_back(var);
    }
    
    // Now set up muParser with all discovered variables
    setupMuParser();
    
    // Test parsing
    m_parser.SetExpr(m_expression);
    m_parser.Eval(); // This will throw if there are syntax errors
    
  } catch (mu::Parser::exception_type &e) {
    throw std::invalid_argument("Expression parsing error: " + std::string(e.GetMsg()));
  } catch (const std::exception& e) {
    throw std::invalid_argument("Expression parsing error: " + std::string(e.what()));
  }
}

void Expression::setupMuParser() const {
  if (m_parser_initialized) {
    return; // Already initialized
  }

  // Clear only variables, keep built-in functions and operators
  m_parser.ClearVar();
  
  // Define the independent variable
  m_parser.DefineVar("x", &m_x_value);
  
  // Resize parameter values vector if needed
  if (m_param_values.size() != m_param_names.size()) {
    m_param_values.resize(m_param_names.size(), 0.0);
  }
  
  // Define parameter variables
  for (size_t i = 0; i < m_param_names.size(); ++i) {
    m_parser.DefineVar(m_param_names[i], &m_param_values[i]);
  }
  
  // Set the expression
  m_parser.SetExpr(m_expression);
  m_parser_initialized = true;
}

double Expression::evaluate(const std::vector<double>& params, double Ecm_over_mref) const {
  try {
    // Set up muParser if not already done
    setupMuParser();
    
    // Set the independent variable directly to Ecm_over_mref
    m_x_value = Ecm_over_mref;
    
    // Set parameter values
    for (size_t i = 0; i < m_param_names.size(); ++i) {
      if (i < m_param_indices.size() && m_param_indices[i] < params.size()) {
        if (i < m_param_values.size()) {
          m_param_values[i] = params[m_param_indices[i]];
        }
      } else {
        if (i < m_param_values.size()) {
          m_param_values[i] = 0.0; // Set to zero as fallback
        }
      }
    }
    
    // Evaluate the expression using muParser
    double result = m_parser.Eval();
    
    return result;
    
  } catch (mu::Parser::exception_type &e) {
    return std::numeric_limits<double>::quiet_NaN();
  } catch (const std::exception& e) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

void Expression::output(XMLHandler& xmlout) const {
  xmlout.set_root("Expression");
  xmlout.put_child("String", m_expression);
}

void Expression::Kinitialize(const KElementInfo& kelem,
                                         std::map<KFitParamInfo, uint>& paramindices) {
  m_param_indices.resize(m_param_names.size());
  
  for (size_t i = 0; i < m_param_names.size(); ++i) {
    // Create a KFitParamInfo for this parameter using the string constructor
    KFitParamInfo kpinfo(m_param_names[i], kelem);
    
    auto it = paramindices.find(kpinfo);
    if (it != paramindices.end()) {
      m_param_indices[i] = it->second;
    } else {
      m_param_indices[i] = paramindices.size();
      paramindices.insert(std::make_pair(kpinfo, m_param_indices[i]));
    }
  }
}

// ***************************************************************************************
