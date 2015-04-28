#ifndef PHG4BlockDetector_h
#define PHG4BlockDetector_h

#include "g4main/PHG4Detector.h"

#include <Geant4/globals.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>


class G4Material;
class G4Box;
class G4LogicalVolume;
class G4Region;
class G4VPhysicalVolume;

class PHG4BlockDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4BlockDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK", const int lyr = 0 );

  //! destructor
  virtual ~PHG4BlockDetector( void )
  {}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  bool IsInBlockActive(G4VPhysicalVolume*) const;
  bool IsInBlock(G4VPhysicalVolume*) const;
  //@}

  void SetSize(const G4double sizex, const G4double sizey, const G4double sizez)
  {dimension[0] = sizex*cm; dimension[1] = sizey*cm; dimension[2] = sizez*cm;}

  void SetMaterial(const std::string &mat) {material = mat;}
  void SetCenterZ(const G4double center_z) {center_in_z = center_z*cm;}
  void SetCenter(const G4double center_x, const G4double center_y, const G4double center_z)
  {
    center_in_x = center_x*cm;
    center_in_y = center_y*cm;
    center_in_z = center_z*cm;
  }

  void SetZRot(const G4double z_angle) {z_rot = z_angle*rad;}

  virtual G4UserSteppingAction* GetSteppingAction() 
  { 
    if ( _region )
      return _region->GetRegionalSteppingAction();
    else return 0;
  }

  void SetActive(const int i = 1) {active = i;}
  int IsActive() const {return active;}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  void BlackHole(const int i=1) {blackhole = i;}
  int IsBlackHole() const {return blackhole;}

  private:

 
  G4Material* TrackerMaterial;
  G4double TrackerThickness;

  G4Material* InactiveMaterial;
  G4double InactiveThickness;

  G4Box* block_solid;
  G4LogicalVolume* block_logic;
  G4VPhysicalVolume* block_physi;


  G4double dimension[3];
  G4String material;
  G4int num_planes_;
  G4double center_in_x;
  G4double center_in_y;
  G4double center_in_z;
  G4double z_rot;

  G4Region* _region;
  int active;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
  
};

#endif
