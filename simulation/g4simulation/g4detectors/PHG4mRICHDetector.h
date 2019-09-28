// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*===============================================================*
 *                        March 2nd 2017                         *
 *         mRICH Detector created by Cheuk-Ping Wong @GSU        *
 *===============================================================*/
#ifndef G4DETECTORS_PHG4MRICHDETECTOR_H
#define G4DETECTORS_PHG4MRICHDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Colour.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Types.hh>        // for G4double, G4int

#include <map>                      // for map
#include <string>

class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;
class PHParameters;
class PHCompositeNode;
class PHG4Subsystem;

//___________________________________________________________________________
class PHG4mRICHDetector: public PHG4Detector
{

 public:
  
  //! constructor
  PHG4mRICHDetector(PHG4Subsystem* subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr = 0);
  
  //! destructor
  virtual ~PHG4mRICHDetector();
  
  //! construct
  virtual void ConstructMe( G4LogicalVolume* world );
  
  //name volume accessors
  //bool IsInBlock(G4VPhysicalVolume*) const;
  int IsInmRICH(G4VPhysicalVolume*) const;

  //void BlackHole(const int i=1) {blackhole = i;}
  //int IsBlackHole() const {return blackhole;}

  void SetActive(const int i = 1)
  {
    active = i;
  }

  void SetAbsorberActive(const int i = 1)
  {
    absorberactive = i;
  }

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

  enum
  {
    SENSOR = 1,
    AEROGEL = 0,
    INACTIVE = -100
  };

  enum DetectorSetUp
  {
    kSingle_Modular = -1,
    kHSector_EWall = 0,
    kHSector = 1,
    kEWall = 2,
    kHWall = 3,
    kHWall_EWall = 4
  };

 private:
  class mRichParameter;
  class BoxPar;
  class PolyPar;
  class LensPar;

  PHParameters *params;
  
  G4VPhysicalVolume* build_box(BoxPar* par, G4LogicalVolume* motherLV);
  G4VPhysicalVolume* build_polyhedra(PolyPar* par, G4LogicalVolume* motherLV);

  G4LogicalVolume* Construct_a_mRICH(G4LogicalVolume* logicWorld);//, int detectorSetup);    //single mRICH
  G4VPhysicalVolume* build_holderBox(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);
  void build_foamHolder(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);
  void build_aerogel(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV);
  void build_lens(LensPar* par, G4LogicalVolume* motherLV);
  void build_mirror(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV);
  void build_sensor(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);

  void build_mRICH_wall_hside(G4LogicalVolume* space);
  void build_mRICH_wall_eside(G4LogicalVolume* space);
  void build_mRICH_sector(G4LogicalVolume* logicWorld, int numSector);
  
  int layer;
  int active;
  int absorberactive;
  //int blackhole;
  std::string superdetector;
  G4VPhysicalVolume *mRICH_PV;         //physical volume of detector box of single module  
  G4VPhysicalVolume *sensor_PV[4];     //physical volume of sensors the sensitive components

  std::map<const G4VPhysicalVolume*, int> sensor_vol; // physical volume of senseors
  std::map<const G4VPhysicalVolume*, int> aerogel_vol; // physical volume of senseors
};
//___________________________________________________________________________
class PHG4mRICHDetector::mRichParameter
{
 private:
  BoxPar* holderBox;
  BoxPar* hollowVolume;
  BoxPar* foamHolderBox;
  PolyPar* foamHolderPoly;
  BoxPar* aerogel;
  LensPar* fresnelLens;
  PolyPar* mirror;
  BoxPar* glassWindow;
  BoxPar* sensor;
  PolyPar* readout;

 public:
  mRichParameter();
  ~mRichParameter();

  void SetPar_glassWindow(int i, G4double x, G4double y);
  void SetPar_sensor(int i, G4double x, G4double y);
  BoxPar* GetBoxPar(std::string componentName);
  LensPar* GetLensPar(std::string componentName);
  PolyPar* GetPolyPar(std::string componentName);

};
//___________________________________________________________________________
class PHG4mRICHDetector::BoxPar
{
 public:
  std::string name;
  G4double halfXYZ[3];
  G4ThreeVector pos;
  G4Material* material;
  int sensitivity;

  G4Colour color;
  bool visibility;
  bool wireframe;
  bool surface;

  BoxPar();
  ~BoxPar();
};
//___________________________________________________________________________
class PHG4mRICHDetector::PolyPar
{
 public:
  std::string name;
  G4ThreeVector pos;
  G4double start;
  G4double theta;
  G4int numSide;
  G4int num_zLayer;
  G4double z[4];                      //max. layer is 4                                                                                                     
  G4double rinner[4];
  G4double router[4];
  G4Material* material;
  int sensitivity;

  G4Colour color;
  bool visibility;
  bool wireframe;
  bool surface;

  PolyPar();
  ~PolyPar(){}
};
//___________________________________________________________________________
class PHG4mRICHDetector::LensPar
{
 public:
  std::string name;
  G4double n;
  G4double f;
  G4double diameter;
  G4double eff_diameter;
  G4double centerThickness;
  G4double grooveWidth;

  G4double halfXYZ[3];
  G4ThreeVector pos;
  G4Material* material;
  int sensitivity;

  G4Colour color;
  bool visibility;
  bool wireframe;
  bool surface;
  
  LensPar();
  ~LensPar(){}

  void Set_halfXYZ(G4double halfX,G4double grooveDensity);
  G4double GetSagita(G4double r);

};
//___________________________________________________________________________
#endif
