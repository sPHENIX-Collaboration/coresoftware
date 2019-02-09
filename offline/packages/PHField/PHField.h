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
  {
  }
  virtual ~PHField() {}
  //! access field value
  //! Follow the convention of G4ElectroMagneticField
  //! @param[in]  Point   space time coordinate. x, y, z, t in Geant4/CLHEP units
  //! @param[out] Bfield  field value. In the case of magnetic field, the order is Bx, By, Bz in in Geant4/CLHEP units
  virtual void GetFieldValue(
      const double Point[4],
      double *Bfield) const = 0;

  void Verbosity(const int i) { m_Verbosity = i; }
  int Verbosity() const { return m_Verbosity; }

 protected:
  int m_Verbosity;
};

#endif
