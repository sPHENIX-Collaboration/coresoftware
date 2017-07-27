/*===============================================================*
 *                        March 2nd 2017                         *
 *         mRICH Detector created by Cheuk-Ping Wong @GSU        *
 *===============================================================*/
#ifndef PHG4mRICHDetector_h
#define PHG4mRICHDetector_h

#include <string>

#include <g4main/PHG4Detector.h>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Colour.hh>

class G4LogicalVolume;
class PHG4Parameters;
class G4VPhysicalVolume;
class G4Material;

using namespace std;
//___________________________________________________________________________
class PHG4mRICHDetector: public PHG4Detector
{

 public:
  
  //! constructor
  PHG4mRICHDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam="BLOCK", const int lyr = 0);
  
  //! destructor
  virtual ~PHG4mRICHDetector( void ) {}
  
  //! construct
  virtual void Construct( G4LogicalVolume* world );
  
  //name volume accessors
  //bool IsInBlock(G4VPhysicalVolume*) const;
  bool IsInmRICH(G4VPhysicalVolume*) const;

  //void BlackHole(const int i=1) {blackhole = i;}
  //int IsBlackHole() const {return blackhole;}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  int get_Layer() const {return layer;}

 private:
  class mRichParameter;
  class BoxPar;
  class PolyPar;
  class LensPar;

  PHG4Parameters *params;
  
  G4VPhysicalVolume* build_box(BoxPar* par, G4LogicalVolume* motherLV);
  G4VPhysicalVolume* build_polyhedra(PolyPar* par, G4LogicalVolume* motherLV);

  G4LogicalVolume* build_Space(G4LogicalVolume* logicWorld, G4double (&bowlPar)[4]);

  G4LogicalVolume* Construct_a_mRICH(G4LogicalVolume* logicWorld);    //single mRICH
  G4VPhysicalVolume* build_holderBox(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);
  void build_foamHolder(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);
  void build_aerogel(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV);
  void build_lens(LensPar* par, G4LogicalVolume* motherLV);
  void build_mirror(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV);
  void build_sensor(mRichParameter* detectorParameter,G4LogicalVolume* motherLV);

  void build_mRICH_wall(G4LogicalVolume* space, G4LogicalVolume* a_mRICH, G4double* bowlPar);
  G4double eta2polarAngle(G4double eta);

  int layer;
  //int blackhole;
  std::string superdetector;
  
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

  void SetPar_glassWindow(G4double x, G4double y);
  void SetPar_sensor(G4double x, G4double y);
  BoxPar* GetBoxPar(std::string componentName);
  LensPar* GetLensPar(std::string componentName);
  PolyPar* GetPolyPar(std::string componentName);

};
//___________________________________________________________________________
class PHG4mRICHDetector::BoxPar
{
 public:
  string name;
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
  string name;
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
  ~PolyPar();
};
//___________________________________________________________________________
class PHG4mRICHDetector::LensPar
{
 public:
  G4double n;
  G4double f;
  G4double diameter;
  G4double eff_diameter;
  G4double centerThickness;
  G4double grooveWidth;

  string name;
  G4double halfXYZ[3];
  G4ThreeVector pos;
  G4Material* material;
  int sensitivity;

  G4Colour color;
  bool visibility;
  bool wireframe;
  bool surface;
  
  LensPar();
  ~LensPar();

  void Set_halfXYZ(G4double halfX,G4double grooveDensity);
  G4double GetSagita(G4double r);

};
//___________________________________________________________________________
#endif
