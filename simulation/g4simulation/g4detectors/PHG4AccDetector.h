#ifndef PHG4AccDetector_h
#define PHG4AccDetector_h

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <set>

class G4EllipticalTube;
class G4LogicalVolume;
class G4Material;
class G4Tubs;
class G4VPhysicalVolume;

class PHG4AccDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4AccDetector( PHCompositeNode *Node, const std::string &dnam, const int layer = 0 );

  //! destructor
  virtual ~PHG4AccDetector( void )
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

  void SetTilt(const double tilt)        {_sciTilt  = tilt;}
  void SetScintWidth(const double wid)   {_sciWidth = wid*cm;}
  void SetTungstenWidth(const double wid)   {_tungstenWidth = wid*cm;}
  void SetUndulMag(const double unmag) {_unmag = unmag*cm;}
  void SetNumScint(const int num)        {_sciNum   = num;}

  void SetThickness(const G4double dbl) {TrackerThickness = dbl*cm;}
  void SetMaterial(const std::string &mat) {material = mat;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  void SetDetectorType(const std::string &typ) {detector_type = typ;}
  int IsInCylinderActive(const G4VPhysicalVolume*);
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}
  G4UserSteppingAction* GetSteppingAction() 
  { 
    if ( _region )
      return _region->GetRegionalSteppingAction();
    else return 0;
  }
  void UseDrawing(const int i=1) {usedrawing = i;}

  void Print(const std::string &what = "ALL") const;

  static int INACTIVE;

  private:

  G4Material* TrackerMaterial;
  G4double    TrackerThickness;

  G4EllipticalTube* ell;

  G4Tubs* cylinder_solid;
  G4LogicalVolume* cylinder_logic;
  G4LogicalVolume* scinti_logic;
  G4VPhysicalVolume* cylinder_physi;
  std::map<const G4VPhysicalVolume*, int> scinti_vol;


  G4double radius;
  G4double length;
  G4double xpos,ypos,zpos;
  G4double _sciTilt;
  G4double _sciWidth;
  G4double _tungstenWidth;
  G4double _unmag;
  G4int    _sciNum;
  G4Region* _region;
  std::string material;
  int active;
  int absorberactive;
  int layer;
  std::string detector_type;
  std::string superdetector;
  double undulations;
  int usedrawing;
};

#endif
