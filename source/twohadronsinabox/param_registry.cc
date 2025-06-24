#include "param_registry.h"
#include <charconv>
#include <regex>

using namespace std;

// ****************************************************************
// *                ParameterNameRegistry Implementation         *
// ****************************************************************

ParameterNameRegistry& ParameterNameRegistry::getInstance() {
  static ParameterNameRegistry instance;
  return instance;
}

uint ParameterNameRegistry::registerParameter(const std::string& name) {

  // Check if name is already registered
  auto name_it = m_name_to_hash.find(name);
  if (name_it != m_name_to_hash.end()) {
    return name_it->second;
  }

  // Compute hash for new parameter name
  uint hash = computeHash(name);

  // Handle hash collisions by finding an unused hash
  while (m_hash_to_name.find(hash) != m_hash_to_name.end()) {
    // If collision and the existing name is different, increment hash
    if (m_hash_to_name[hash] != name) {
      hash++;
    } else {
      // Same name already exists with this hash
      return hash;
    }
  }

  // Register the new parameter
  m_hash_to_name[hash] = name;
  m_name_to_hash[name] = hash;

  return hash;
}

std::string ParameterNameRegistry::getParameterName(uint hash) const {
  auto it = m_hash_to_name.find(hash);
  if (it != m_hash_to_name.end()) {
    return it->second;
  }
  return std::string(); // Return empty string if not found
}

bool ParameterNameRegistry::hasHash(uint hash) const {
  return m_hash_to_name.find(hash) != m_hash_to_name.end();
}

void ParameterNameRegistry::clear() {
  m_hash_to_name.clear();
  m_name_to_hash.clear();
}

size_t ParameterNameRegistry::size() const {
  return m_hash_to_name.size();
}

uint ParameterNameRegistry::computeHash(const std::string& name) {
  // Use the same hashing algorithm as in KFitParamInfo for expressions
  uint hash = 0;
  for (char c : name) {
    hash = hash * 31 + static_cast<uint>(c);
  }
  // Mask to 28 bits to match KFitParamInfo storage constraint
  return hash & 0xFFFFFFFu;
}

std::string ParameterNameRegistry::getParameterNameFromMCObsName(const std::string& obs_name) const {
  uint hash = parseParamHashFromMCObsName(obs_name);
  if (hash == 0) {
    return {}; // Not a string expression parameter
  }
  return getParameterName(hash);
}

uint ParameterNameRegistry::getHashFromMCObsName(const std::string& obs_name) const {
  return parseParamHashFromMCObsName(obs_name);
}


uint ParameterNameRegistry::parseParamHashFromMCObsName(std::string_view obs_name) {
  const std::string kPrefix = "KStrExpr(";
  
  // Convert string_view to string for compatibility
  std::string obs_str(obs_name);
  
  if (obs_str.size() < kPrefix.size() || 
      obs_str.substr(0, kPrefix.size()) != kPrefix) {
    return 0; // Not a valid KStrExpr name
  }

  // Position right after the '('
  std::string rest = obs_str.substr(kPrefix.size());
  
  // Find the first comma or closing parenthesis
  size_t end_pos = rest.find_first_of(",)");
  if (end_pos == std::string::npos) return 0; // Invalid format

  // Extract just the hash part
  std::string hash_part = rest.substr(0, end_pos);
  if (hash_part.empty()) return 0; // Empty hash

  // Convert the number up to the first ',' or ')' using more compatible approach
  try {
    size_t pos = 0;
    uint hash = std::stoul(hash_part, &pos);
    
    // Check if entire string was consumed
    if (pos != hash_part.size()) return 0;
    
    return hash;
  } catch (const std::exception&) {
    return 0; // parse failure
  }
}