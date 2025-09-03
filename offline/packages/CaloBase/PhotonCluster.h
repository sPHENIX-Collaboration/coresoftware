#ifndef CALOBASE_PHOTONCLUSTER_H
#define CALOBASE_PHOTONCLUSTER_H

#include <phool/phool.h>
#include <Rtypes.h>  // for ROOT dictionary macro definitions
#include <iostream>
#include <limits>

//! Interface mixin providing photon-specific augmentation to a calorimeter cluster.
//!
//! Authors:
//!  - S. Li <sli7@bnl.gov>
//!
//! Notes:
//!  - Energy, position and isolation ET live in the underlying RawCluster implementation.
//!  - This interface only covers photon ID related quantities (conversion + shower shape) and
//!    convenience hooks for ID logic/printing. (conversion flag is a placeholder)
//!  - Keep it lean: no data members here – storage resides in concrete subclasses.
class PhotonCluster 
{
 public:
  //! Virtual destructor for proper cleanup via base pointer
  virtual ~PhotonCluster();

  //! @name Photon Property Getters
  //! @{
  //! Implement in derived classes. Default versions only warn and return sentinel values.
  virtual float get_conversion_probability() const
  {
    PHOOL_VIRTUAL_WARN("get_conversion_probability()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual bool is_converted() const
  {
    PHOOL_VIRTUAL_WARN("is_converted()");
    return false;
  }

  virtual float get_shower_shape_parameter(const std::string& /*name*/) const
  {
    PHOOL_VIRTUAL_WARN("get_shower_shape_parameter()");
    return std::numeric_limits<float>::quiet_NaN();
  }
  //! @}

  //! @name Photon Analysis Methods  
  //! @{
  virtual bool pass_photon_cuts() const
  {
    PHOOL_VIRTUAL_WARN("pass_photon_cuts()");
    return false;
  }

  virtual void identify_photon(std::ostream& /*os*/ = std::cout) const
  {
    PHOOL_VIRTUAL_WARN("identify_photon()");
  }
  //! @}

  //! @name Photon Property Setters
  //! @{
  virtual void set_conversion_probability(const float /*prob*/)
  {
    PHOOL_VIRTUAL_WARN("set_conversion_probability()");
  }

  virtual void set_converted(const bool /*converted*/)
  {
    PHOOL_VIRTUAL_WARN("set_converted()");
  }

  virtual void set_shower_shape_parameter(const std::string& /*name*/,  const float /*shape*/)
  {
    PHOOL_VIRTUAL_WARN("set_shower_shape_parameter()");
  }
  //! @}

  //! @name Lifecycle / validation helpers
  //! @{
  virtual void reset_photon_properties()
  {
    PHOOL_VIRTUAL_WARN("reset_photon_properties()");
  }
  
  virtual bool is_valid_photon() const
  {
    PHOOL_VIRTUAL_WARN("is_valid_photon()");
    return false;
  }
  //! @}

 protected:
  // Prevent direct instantiation; allow construction by derived classes only.
  PhotonCluster() = default;

  // Interface-only class (no PHObject inheritance) -> no overrides of ROOT virtuals
  // Use plain ClassDef to generate dictionary if needed for containers holding PhotonCluster*.
  ClassDef(PhotonCluster, 1)
};

#endif
