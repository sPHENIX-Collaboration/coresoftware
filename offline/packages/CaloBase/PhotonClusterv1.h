#ifndef CALOBASE_PHOTONCLUSTERV1_H
#define CALOBASE_PHOTONCLUSTERV1_H

#include "PhotonCluster.h"
#include "RawClusterv1.h"

#include <map>
#include <string>

//! PhotonClusterv1 - derives from both PhotonCluster and RawClusterv1
//! @warning Multiple inheritance used - be careful with method resolution
class PhotonClusterv1 : public PhotonCluster, public RawClusterv1
{
 public:
  PhotonClusterv1() = default;
  
  //! @warning Virtual destructor override - ensures proper cleanup
  ~PhotonClusterv1() override = default;

  //! Copy constructor from existing RawClusterv1 object
  //! @warning This will copy all RawClusterv1 data but initialize photon properties to defaults
  explicit PhotonClusterv1(const RawClusterv1& rawcluster);

  //! Copy constructor 
  PhotonClusterv1(const PhotonClusterv1& other) = default;
  
  //! Assignment operator
  PhotonClusterv1& operator=(const PhotonClusterv1& other) = default;

  //! @name PHObject Interface Overrides
  //! @{
  //! @warning Override methods from RawClusterv1 - virtual dispatch applies
  void Reset() override;
  PHObject* CloneMe() const override { return new PhotonClusterv1(*this); }
  int isValid() const override;
  void identify(std::ostream& os = std::cout) const override;
  //! @}

  //! @name PhotonCluster Virtual Method Implementations
  //! @{
  //! @warning These methods override virtual functions from PhotonCluster base
  float get_conversion_probability() const override { return m_conversion_prob; }
  bool is_converted() const override { return m_is_converted; }
  // RawClusterv1 already provides energy and isolation energy accessors
  // No single default shower shape parameter anymore
  bool pass_photon_cuts() const override;
  void identify_photon(std::ostream& os = std::cout) const override;
  bool is_valid_photon() const override;
  void reset_photon_properties() override;
  //! @}

  //! @name PhotonCluster Setter Implementations
  //! @{
  //! @warning Virtual override - can be further overridden in derived classes
  void set_conversion_probability(const float prob) override { m_conversion_prob = prob; }
  void set_converted(const bool converted) override { m_is_converted = converted; }
  // Extra helper to set named shower shape parameters
  void set_shower_shape_parameter(const std::string& name, float value) override { m_shower_shapes[name] = value; }
  float get_shower_shape_parameter(const std::string& name) const override;
  const std::map<std::string,float>& get_all_shower_shapes() const { return m_shower_shapes; }
  //! @}

 private:
  //! @warning Photon-specific data members - memory managed only in this derived class
  //! @warning Base class PhotonCluster has NO data members by design
  // Photon energy and isolation energy now sourced from RawClusterv1
  float m_conversion_prob{0.0f};        //!< Probability of photon conversion
  bool m_is_converted{false};           //!< Conversion flag
  std::map<std::string,float> m_shower_shapes; //!< Named shower shape parameters

  ClassDefOverride(PhotonClusterv1, 1)  //!< ROOT dictionary generation
};

#endif
