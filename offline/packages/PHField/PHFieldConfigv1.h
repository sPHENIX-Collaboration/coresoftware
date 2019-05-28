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

  virtual ~PHFieldConfigv1() {}

  /// Virtual copy constructor.
  virtual PHObject*
  clone() const;

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const;

  /// Clear Content
  virtual void
  Reset() {}

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

  ClassDef(PHFieldConfigv1, 3)
};

#endif
