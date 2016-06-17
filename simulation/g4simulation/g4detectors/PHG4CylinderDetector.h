#ifndef PHG4CylinderDetector_h
#define PHG4CylinderDetector_h

#include "g4main/PHG4Detector.h"

#include <Geant4/globals.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <float.h>

class G4Material;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

class PHG4CylinderDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4CylinderDetector( PHCompositeNode *Node, const std::string &dnam, const int layer = 0 );

  //! destructor
  virtual ~PHG4CylinderDetector( void )
  {}

  //! construct
  void Construct( G4LogicalVolume* world );

  void SetRadius(const G4double dbl) {radius = dbl*cm;}
  void SetLength(const G4double dbl) {length = dbl*cm;}
  void SetPosition(const G4double x, const G4double y, const G4double z)
  {
    xpos = x*cm;
    ypos = y*cm;
    zpos = z*cm;
  }
  void SetThickness(const G4double dbl) {TrackerThickness = dbl*cm;}
  void SetMaterial(const std::string &mat) {material = mat;}
  void SetActive(const int i = 1) {active = i;}
  int IsActive() const {return active;}
  void SetReducedTruthInfo(const bool b) {reduced = b;}
  bool GetReducedTruthInfo() {return reduced;}
  void SetDetectorType(const std::string &typ) {detector_type = typ;}
  bool IsInCylinderActive(const G4VPhysicalVolume*) const;
  bool IsInCylinder(const G4VPhysicalVolume*) const;
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}
  G4UserSteppingAction* GetSteppingAction() 
  { 
    if ( _region )
      return _region->GetRegionalSteppingAction();
    else return 0;
  }

  void BlackHole(const int i=1,
		 const double tmin = DBL_MAX,
		 const double tmax = -1.0*DBL_MAX) {
    blackhole = i; blackhole_tmin = tmin; blackhole_tmax = tmax;
  }
  int IsBlackHole() const {return blackhole;}
  double GetBlackHoleTMin() {return blackhole_tmin;}
  double GetBlackHoleTMax() {return blackhole_tmax;}
  
  private:

  G4Material* TrackerMaterial;
  G4double TrackerThickness;

  G4Tubs* cylinder_solid;
  G4LogicalVolume* cylinder_logic;
  G4VPhysicalVolume* cylinder_physi;


  G4double radius;
  G4double length;
  G4double xpos,ypos,zpos;
  G4Region* _region;
  std::string material;
  int active;
  bool reduced;
  int layer;
  int blackhole;
  double blackhole_tmin;
  double blackhole_tmax;
  std::string detector_type;
  std::string superdetector;
};

#endif
