#ifndef PARAM_REGISTRY_H
#define PARAM_REGISTRY_H

#include <string>
#include <algorithm>


// ***********************************************************************
// *                                                                     *
// *   Registry class for storing parameter names by their hashes        *
// *                                                                     *
// ***********************************************************************

//  The "ParameterNameRegistry" class provides a centralized registry for
//  storing and retrieving parameter names by their hash values. This avoids
//  storing duplicate parameter name strings in multiple objects while still
//  allowing retrieval of the original names when needed.


class ParameterNameRegistry {
public:
  // Get the singleton instance
  static ParameterNameRegistry& getInstance();

  // Register a parameter name and return its hash
  // If the name is already registered, returns the existing hash
  uint registerParameter(const std::string& name);

  std::string getParameterName(uint hash) const;

  // Check if a hash exists in the registry
  bool hasHash(uint hash) const;

  // Clear the registry (mainly for testing purposes)
  void clear();

  // Get the number of registered parameters
  size_t size() const;

  // Helper method to get parameter name directly from KFitParamInfo object
  // Returns empty string if not a string expression parameter or hash not found
  std::string getParameterNameFromMCObsName(const std::string& obs_name) const;

  uint getHashFromMCObsName(const std::string& obs_name) const;

private:
  ParameterNameRegistry() = default;
  ~ParameterNameRegistry() = default;

  // Disable copy/move constructors and assignment operators
  ParameterNameRegistry(const ParameterNameRegistry&) = delete;
  ParameterNameRegistry& operator=(const ParameterNameRegistry&) = delete;
  ParameterNameRegistry(ParameterNameRegistry&&) = delete;
  ParameterNameRegistry& operator=(ParameterNameRegistry&&) = delete;

  // Compute hash of a parameter name
  static uint computeHash(const std::string& name);

  static uint parseParamHashFromMCObsName(std::string_view obs_name);

  std::unordered_map<uint, std::string> m_hash_to_name;
  std::unordered_map<std::string, uint> m_name_to_hash;
};

#endif //PARAM_REGISTRY_H
