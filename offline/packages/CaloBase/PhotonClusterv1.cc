#include "PhotonClusterv1.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

PhotonClusterv1::PhotonClusterv1(const RawCluster& rc)
  : RawClusterv1([&rc]() -> const RawClusterv1& {
      if (const auto* rc1 = dynamic_cast<const RawClusterv1*>(&rc))
      {
        return *rc1;
      }
      throw std::runtime_error("PhotonClusterv1 requires RawClusterv1 (or derived).");
    }())
{
  if (const auto* photon = dynamic_cast<const PhotonClusterv1*>(&rc))
  {
    m_shower_shapes = photon->get_all_shower_shapes();
  }
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
  m_shower_shapes.clear();
  return;
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
  // List all named shower shapes
  for (const auto& kv : m_shower_shapes)
  {
    os << "    shape[" << kv.first << "]: " << kv.second << std::endl;
  }
  os << "------------------------------------------" << std::endl;
}

bool PhotonClusterv1::pass_photon_cuts() const
{
  // @warning: These are example cuts - customize based on your analysis needs
  std::cout << "PhotonClusterv1::pass_photon_cuts() is currently unimplemented" << std::endl;
  // Minimum energy cut

  // Shower shape cut: if a shape named "core" exists, apply cut
  // auto it_core = m_shower_shapes.find("core");
  // if (it_core != m_shower_shapes.end()) {
  //  if (it_core->second > 0.3f) return false;
  //}


  // @warning: Add more sophisticated photon ID cuts as needed
  // Consider using cluster properties like get_ecore(), get_prob(), etc.

  return true;
}

float PhotonClusterv1::get_shower_shape_parameter(const std::string& name) const
{
  auto it = m_shower_shapes.find(name);
  if (it != m_shower_shapes.end())
  {
    return it->second;
  }
  return std::numeric_limits<float>::quiet_NaN();
}