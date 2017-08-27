// $Id: $

/*!
 * \file PHFieldConfig_v1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFieldConfig_v1_H_
#define PHFieldConfig_v1_H_

#include <vector>
#include "PHFieldConfig.h"

class TGeoVolume;
class TGeoManager;

/*!
 * \brief PHFieldConfig_v1 impliments field configuration information */
class PHFieldConfig_v1 : public PHFieldConfig
{
 public:
  PHFieldConfig_v1(
      FieldConfigTypes field_config,
      std::string& filename,
      double magfield_rescale = 1.);
  virtual ~PHFieldConfig_v1();

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void
  Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int
  isValid() const;

  FieldConfigTypes get_field_config() const
  {
    return field_config_;
  }
  void set_field_config(FieldConfigTypes fieldConfig)
  {
    field_config_ = fieldConfig;
  }

  const std::string& get_filename() const
  {
    return filename_;
  }

  void set_filename(const std::string& filename)
  {
    filename_ = filename;
  }

  double get_magfield_rescale() const
  {
    return magfield_rescale_;
  }

  void set_magfield_rescale(double magfieldRescale)
  {
    magfield_rescale_ = magfieldRescale;
  }

 protected:
  FieldConfigTypes field_config_;
  std::string filename_;
  double magfield_rescale_;

  ClassDef(PHFieldConfig_v1, 3)
};

#endif /* PHFieldConfig_v1_H_ */
