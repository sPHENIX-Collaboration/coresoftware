// $Id: $

/*!
 * \file PHFieldConfig_v1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFIELD_PHFIELDCONFIGV1_H
#define PHFIELD_PHFIELDCONFIGV1_H

#include "PHFieldConfig.h"

#include <iostream>
#include <string>

class PHObject;

/*!
 * \brief PHFieldConfigv1 implements field configuration information for input a field map file */
class PHFieldConfigv1 : public PHFieldConfig
{
 public:
  PHFieldConfigv1(
      FieldConfigTypes field_config,
      const std::string& filename,
      double magfield_rescale = 1.);

  //! default constructor for ROOT file IO
  PHFieldConfigv1()
    : PHFieldConfigv1(kFieldInvalid, "INVALID FILE")
  {
  }

  ~PHFieldConfigv1() override {}

  /// Virtual copy constructor.
  PHObject* CloneMe() const override { return new PHFieldConfigv1(*this); }

  /** identify Function from PHObject
   @param os Output Stream
   */
  void
  identify(std::ostream& os = std::cout) const override;

  /// Clear Content
  void Reset() override {}

  /// isValid returns non zero if object contains vailid data
  int
  isValid() const override;

  FieldConfigTypes get_field_config() const override
  {
    return field_config_;
  }
  void set_field_config(FieldConfigTypes fieldConfig) override
  {
    field_config_ = fieldConfig;
  }

  const std::string& get_filename() const override
  {
    return filename_;
  }

  void set_filename(const std::string& filename) override
  {
    filename_ = filename;
  }

  double get_magfield_rescale() const override
  {
    return magfield_rescale_;
  }

  void set_magfield_rescale(double magfieldRescale) override
  {
    magfield_rescale_ = magfieldRescale;
  }

 protected:
  FieldConfigTypes field_config_;
  std::string filename_;
  double magfield_rescale_;

  ClassDefOverride(PHFieldConfigv1, 3)
};

#endif
