#ifndef PHFIELD_PHFIELD_H
#define PHFIELD_PHFIELD_H

// units of this class. To convert internal value to Geant4/CLHEP units for fast access

//! \brief transient object for field storage and access
class PHField
{
 public:
  //! constructor
  explicit PHField(const int verb = 0)
    : m_Verbosity(verb)
  {}

  //! destructor
  virtual ~PHField() = default;

  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  virtual void GetFieldValue(
      const double Point[4],
      double *Bfield) const = 0;

  //! un-cached version of field accessor.
  /* used for multi-threading. By default, the same as GetFieldValue */
  virtual void GetFieldValue_nocache(
      const double Point[4],
      double *Bfield) const
  { return GetFieldValue( Point, Bfield ); }

  //! verbosity
  void Verbosity(const int i) { m_Verbosity = i; }

  //! verbosity
  int Verbosity() const { return m_Verbosity; }

 protected:

  //! verbosity
  int m_Verbosity = 0;
};

#endif
