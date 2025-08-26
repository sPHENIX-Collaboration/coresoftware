#include "PhotonClusterv1.h"
#include <iostream>
#include <limits>
#include <map>
#include <string>

PhotonClusterv1::PhotonClusterv1(const RawClusterv1& rawcluster)
  : PhotonCluster()
  , RawClusterv1(rawcluster)
  , m_conversion_prob(0.0f)
  , m_is_converted(false)
  , m_shower_shapes()
{
  // No default shower shape inserted; user/algorithm must set if needed
  // Photon energy now same as cluster energy via get_energy()
  // Isolation energy sourced from RawClusterv1 get_et_iso()
}

void PhotonClusterv1::Reset()
{
  // @warning: Call base class reset first to avoid issues with virtual dispatch
  RawClusterv1::Reset();
  
  // Reset photon-specific members
  reset_photon_properties();
}

void PhotonClusterv1::reset_photon_properties()
{
  // Reset photon-related values (cluster energy left to RawClusterv1 Reset caller if needed)
  m_conversion_prob = 0.0f;
  m_is_converted = false;
  m_shower_shapes.clear();
}

int PhotonClusterv1::isValid() const
{
  // @warning: Multiple inheritance - both base class validations checked
  return RawClusterv1::isValid() && is_valid_photon();
}

bool PhotonClusterv1::is_valid_photon() const
{
  // Valid photon if has positive energy
  return (get_energy() > 0.0f);
}

void PhotonClusterv1::identify(std::ostream& os) const
{
  // @warning: Call base class identify first to maintain output order
  RawClusterv1::identify(os);
  
  // Add photon-specific information
  identify_photon(os);
}

void PhotonClusterv1::identify_photon(std::ostream& os) const
{
  os << "--- PhotonClusterv1 Photon Properties ---" << std::endl;
  os << "  Photon Energy: " << get_energy() << " GeV" << std::endl;
  os << "  Conversion Probability: " << get_conversion_probability() << std::endl;
  os << "  Is Converted: " << (is_converted() ? "Yes" : "No") << std::endl;
  os << "  Isolation Energy: " << get_et_iso() << " GeV" << std::endl;
  // List all named shower shapes
  for (const auto& kv : m_shower_shapes)
  {
    os << "    shape[" << kv.first << "]: " << kv.second << std::endl;
  }
  os << "  Passes Photon Cuts: " << (pass_photon_cuts() ? "Yes" : "No") << std::endl;
  os << "------------------------------------------" << std::endl;
}

bool PhotonClusterv1::pass_photon_cuts() const
{
  // @warning: These are example cuts - customize based on your analysis needs
  std::cout<<"this is currently unimplemented"<<std::endl;
  // Minimum energy cut
  if (get_energy() < 0.5f) {
    return false;  
  }
  
  // Shower shape cut: if a shape named "core" exists, apply cut
  //auto it_core = m_shower_shapes.find("core");
  //if (it_core != m_shower_shapes.end()) {
  //  if (it_core->second > 0.3f) return false;
  //}
  
  // Isolation cut (photons should be isolated)
  if (get_et_iso() > 2.0f) {
    return false;  
  }
  
  // @warning: Add more sophisticated photon ID cuts as needed
  // Consider using cluster properties like get_ecore(), get_prob(), etc.
  
  return true;
}

float PhotonClusterv1::get_shower_shape_parameter(const std::string& name) const
{
  auto it = m_shower_shapes.find(name);
  if (it != m_shower_shapes.end()) return it->second;
  return std::numeric_limits<float>::signaling_NaN();
}