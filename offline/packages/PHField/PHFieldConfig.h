// $Id: $

/*!
 * \file PHFieldConfig.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHFieldConfig_H_
#define PHFieldConfig_H_

#include <phool/PHObject.h>
#include <limits>
#include <string>

/*!
 * \brief PHFieldConfig store field configuration information */
class PHFieldConfig : public PHObject
{
 public:
  virtual ~PHFieldConfig();

  /// Virtual copy constructor.
  virtual PHObject*
  clone() const;

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

  enum FieldConfigTypes
  {
    //! consttant field
    kFieldConstant = 0,
    //! 2D field map expressed in cylindrical coordinates
    kField2D = 2,
    //! 3D field map expressed in cylindrical coordinates
    kField3DCylindrical = 3,
    //! 3D field map expressed in Cartesian coordinates
    Field3DCartesian = 1,

    //! invalid value
    kFieldInvalid = std::numeric_limits<int>::max()
  };

  virtual FieldConfigTypes get_field_config() const
  {
    return kFieldInvalid;
  }

  virtual void set_field_config(FieldConfigTypes fieldConfig)
  {
  }

  virtual const std::string& get_filename() const
  {
    return "INVALID FILE";
  }

  virtual void set_filename(const std::string& filename)
  {
  }

  virtual double get_magfield_rescale() const
  {
    return 1.0;
  }

  virtual void set_magfield_rescale(double magfieldRescale)
  {
  }

 protected:
  //! pure virtual interface class. not for direct use
  PHFieldConfig();

  ClassDef(PHFieldConfig, 1)
};

#endif /* PHFieldConfig_H_ */
