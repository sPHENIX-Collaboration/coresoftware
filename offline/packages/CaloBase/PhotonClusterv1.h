#ifndef CALOBASE_PHOTONCLUSTERV1_H
#define CALOBASE_PHOTONCLUSTERV1_H

#include "RawClusterv1.h"

#include <map>
#include <string>

//! PhotonClusterv1 - extends RawClusterv1 with photon-ID helpers
class PhotonClusterv1 : public RawClusterv1
{
 public:
  PhotonClusterv1() = default;

  ~PhotonClusterv1() override = default;


  explicit PhotonClusterv1(const RawCluster & rc);

  //! Copy constructor
  PhotonClusterv1(const PhotonClusterv1& other) = default;

  //! Assignment operator
  PhotonClusterv1& operator=(const PhotonClusterv1& other) = default;

  void Reset() override;
  PHObject* CloneMe() const override { return new PhotonClusterv1(*this); }
  void identify(std::ostream& os = std::cout) const override;
  bool pass_photon_cuts() const override;
  void identify_photon(std::ostream& os = std::cout) const override;


  void reset_photon_properties() override;

  //! @name PhotonCluster Setter Implementations
  //! @{
  //! @warning Virtual override - can be further overridden in derived classes

  // Extra helper to set named shower shape parameters
  void set_shower_shape_parameter(const std::string& name, float value) override { m_shower_shapes[name] = value; }
  float get_shower_shape_parameter(const std::string& name) const override;
  const std::map<std::string, float>& get_all_shower_shapes() const override { return m_shower_shapes; }
  //! @}

 private:
  //! @warning Photon-specific data members - memory managed only in this derived class
  // Photon energy and isolation energy now sourced from RawCluster
  //float m_conversion_prob{0.0f};                 //!< Probability of photon conversion
  //bool m_is_converted{false};                    //!< Conversion flag
  std::map<std::string, float> m_shower_shapes;  //!< Named shower shape parameters

  ClassDefOverride(PhotonClusterv1, 1)  //!< ROOT dictionary generation
};

#endif
