#ifndef __PHG4FIELDSPHENIX_H__
#define __PHG4FIELDSPHENIX_H__

#include <Geant4/G4MagneticField.hh>


#include <map> 
#include <set>

class PHG4FieldsPHENIX : public G4MagneticField
{
  
 public:
  
  PHG4FieldsPHENIX(const std::string  &fname);
  virtual ~PHG4FieldsPHENIX();
  
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  
 protected:
  std::string filename; 
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  double xstepsize;
  double ystepsize;
  double zstepsize;
  // these are updated in a const method
  // to cache previous values
  mutable double xyz[2][2][2][3];
  mutable double bf[2][2][2][3];
  mutable double xkey_save;
  mutable double ykey_save;
  mutable double zkey_save;
  mutable int cache_hits;
  mutable int cache_misses;
};


#endif // __PHFIELD3D_H
