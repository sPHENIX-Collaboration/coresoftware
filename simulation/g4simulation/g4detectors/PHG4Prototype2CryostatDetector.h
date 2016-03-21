#ifndef PHG4Prototype2CryostatDetector_h
#define PHG4Prototype2CryostatDetector_h

#include "PHG4Parameters.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

class PHG4Prototype2CryostatDetector: public PHG4Detector
{

  public:

  //! constructor
 PHG4Prototype2CryostatDetector( PHCompositeNode *Node,  PHG4Parameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4Prototype2CryostatDetector(){}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInPrototype2Cryostat(G4VPhysicalVolume*) const;
  //@}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  G4LogicalVolume* ConstructAluPlate(G4LogicalVolume* hcalenvelope, const int n);

  protected:
  void AddGeometryNode();
  int ConstructCryostat(G4LogicalVolume* sandwich);
  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);
  PHG4Parameters *params;
  G4AssemblyVolume *cryostatassembly;
  G4TwoVector alu_plate_corner_upper_left[3];
  G4TwoVector alu_plate_corner_upper_right[3];
  G4TwoVector alu_plate_corner_lower_right[3];
  G4TwoVector alu_plate_corner_lower_left[3];
  double alu_z;
  double volume_alu[3];

  int n_alu_plates;
  int active;
  int absorberactive;

  int layer;
  std::string detector_type;
  std::string superdetector;
  std::string scintilogicnameprefix;
};

#endif
