/*===============================================================*
 *                        March 2nd 2017                         *
 *         mRICH Detector created by Cheuk-Ping Wong @GSU        *
 *===============================================================*
 * Even mRICH is a tiny detector, it has nine logical volumes per*
 * modules. To make the code easy to manage, material definition *
 * and components dimensions are written in different functions  *
 * which are not supposed for frequent modification.             *
 *                                                               * 
 * While Construct() and Construct_a_mRIHC(), are two            *
 * handy and simple functions for user to control detector       *
 * construction.                                                 *
 *---------------------------------------------------------------*
 * to ignore a particular component of a single mRICH, comment   *
 * out the build_xxx() function correspond to the component.     *
 *---------------------------------------------------------------*
 * Materials are defined in gmain/PHG4Reco::DefineMaterials      *
 *===============================================================*/
#include "PHG4mRICHDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>  // for PHG4Detector

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalBorderSurface.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialPropertiesTable.hh>  // for G4MaterialProperties...
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SurfaceProperty.hh>  // for dielectric_dielectric
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4ios.hh>  // for G4cout, G4endl

#include <algorithm>  // for fill, max
#include <cmath>  // for floor, sqrt, acos, asin
#include <cstdio>
#include <iostream>  // for operator<<, basic_os...
#include <iterator>  // for begin, end
#include <set>       // for _Rb_tree_const_iterator
#include <string>

class PHCompositeNode;

using namespace std;
using namespace CLHEP;

//_______________________________________________________________
PHG4mRICHDetector::PHG4mRICHDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , params(parameters)
  ,
  //block_physi(nullptr),
  layer(lyr)
  , active(0)
  , absorberactive(0)
  , mRICH_PV(nullptr)
  , sensor_PV{nullptr, nullptr, nullptr, nullptr}
{
}

PHG4mRICHDetector::~PHG4mRICHDetector()
{
}
//_______________________________________________________________
int PHG4mRICHDetector::IsInmRICH(G4VPhysicalVolume* volume) const
{
  if (active && sensor_vol.find(volume) != sensor_vol.end())
  {
    return SENSOR;
  }

  if (active && aerogel_vol.find(volume) != aerogel_vol.end())
  {
    return AEROGEL;
  }

  return INACTIVE;
}
//______________________________________________________________
void PHG4mRICHDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  int subsystemSetup = params->get_int_param("subsystemSetup");
  // -1: single module
  //  0: h-side sectors and e-side wall
  //  1: h-side sectors
  //  2: e-side wall
  //  3: h-side wall
  //  4: h-side wall and e-side wall

  if (subsystemSetup == DetectorSetUp::kSingle_Modular) Construct_a_mRICH(logicWorld);
  if (subsystemSetup == DetectorSetUp::kHSector_EWall)
  {
    build_mRICH_sector(logicWorld, 8);
    build_mRICH_wall_eside(logicWorld);
  }
  if (subsystemSetup == DetectorSetUp::kHSector) build_mRICH_sector(logicWorld, 8);
  if (subsystemSetup == DetectorSetUp::kEWall) build_mRICH_wall_eside(logicWorld);
  if (subsystemSetup == DetectorSetUp::kHWall) build_mRICH_wall_hside(logicWorld);
  if (subsystemSetup == DetectorSetUp::kHWall_EWall)
  {
    build_mRICH_wall_hside(logicWorld);
    build_mRICH_wall_eside(logicWorld);
  }
}
//_______________________________________________________________
G4LogicalVolume* PHG4mRICHDetector::Construct_a_mRICH(G4LogicalVolume* logicWorld)  //, int detectorSetup )
{
  int detectorSetup = params->get_int_param("detectorSetup");

  mRichParameter* parameters = new mRichParameter();

  //--------------------------- skeleton setup ---------------------------//
  /*holder box and hollow volume*/ G4VPhysicalVolume* hollowVol = build_holderBox(parameters, logicWorld);
  /*aerogel                     */ build_aerogel(parameters, hollowVol);
  /*sensor plane                */ build_sensor(parameters, hollowVol->GetLogicalVolume());

  //-------------------------- for full set up ---------------------------//
  if (detectorSetup)
  {  //for full setup
    /*foam holder for aerogel     */ build_foamHolder(parameters, hollowVol->GetLogicalVolume());
    /*lens                        */ build_lens(parameters->GetLensPar("fresnelLens"), hollowVol->GetLogicalVolume());
    /*mirror                      */ build_mirror(parameters, hollowVol);
    /*readout electronics         */ G4VPhysicalVolume* pol = build_polyhedra(parameters->GetPolyPar("readout"), hollowVol->GetLogicalVolume());
    if (!pol)
    {
      cout << "readout electronics not placed" << endl;
    }
  }

  return hollowVol->GetMotherLogical();  //return detector holder box.
                                         //you have more than 1 daugthers,
                                         //but you can only have one mother.
}

//________________________________________________________________________//
PHG4mRICHDetector::BoxPar::BoxPar()
{
  fill(begin(halfXYZ), end(halfXYZ), (G4double) 0 * mm);
  pos = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);
  material = G4Material::GetMaterial("G4_AIR");
  sensitivity = 0;

  color = G4Colour(0, 0, 0, 0);
  visibility = false;
  wireframe = false;
  surface = false;
}
//________________________________________________________________________//
PHG4mRICHDetector::BoxPar::~BoxPar() { ; }
//________________________________________________________________________//
PHG4mRICHDetector::PolyPar::PolyPar()
  : name("")
  , pos(G4ThreeVector(0 * mm, 0 * mm, 0 * mm))
  , start(0.)
  , theta(0.)
  , numSide(0)
  , num_zLayer(0)
  , material(G4Material::GetMaterial("G4_AIR"))
  , sensitivity(0)
  , color(G4Colour(0, 0, 0, 0))
  , visibility(false)
  , wireframe(false)
  , surface(false)
{
  fill(begin(z), end(z), (G4double) 0 * mm);
  fill(begin(rinner), end(rinner), (G4double) 0 * mm);
  fill(begin(router), end(router), (G4double) 0 * mm);
}
//________________________________________________________________________//
PHG4mRICHDetector::LensPar::LensPar()
  : name("")
  , n(0.)
  , f(0.)
  , diameter(0.)
  , eff_diameter(0.)
  , centerThickness(0.)
  , grooveWidth(0.)
  , pos(G4ThreeVector(0 * mm, 0 * mm, 0 * mm))
  , material(G4Material::GetMaterial("G4_AIR"))
  , sensitivity(0.)
  , color(G4Colour(0, 0, 0, 0))
  , visibility(false)
  , wireframe(false)
  , surface(false)

{
  fill(begin(halfXYZ), end(halfXYZ), (G4double) 0 * mm);
}
//________________________________________________________________________//
void PHG4mRICHDetector::LensPar::Set_halfXYZ(G4double halfX, G4double grooveDensity)
{
  halfXYZ[0] = halfX;
  halfXYZ[1] = halfXYZ[0];

  G4double NumberOfGrooves = floor(grooveDensity * (eff_diameter / 2.0));
  G4double Rmin1 = (NumberOfGrooves - 1) * (grooveWidth);
  G4double Rmax1 = (NumberOfGrooves - 0) * (grooveWidth);
  halfXYZ[2] = (GetSagita(Rmax1) - GetSagita(Rmin1) + centerThickness) / 2.0;
}
//________________________________________________________________________//
G4double PHG4mRICHDetector::LensPar::GetSagita(G4double r)
{
  G4double Conic = -1.0;  // original:
  //G4int lens_type = 3;
  G4int lens_type = 5;
  G4double Curvature;
  G4double Aspher[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  if (lens_type == 1)
  {
    Curvature = 0.00437636761488 / mm;
    Aspher[0] = 4.206739256e-05 / (mm);
    Aspher[1] = 9.6440152e-10 / (mm3);
    Aspher[2] = -1.4884317e-15 / (mm2 * mm3);
  }
  else if (lens_type == 2)
  {
    Curvature = 0.0132 / mm;  // r=77mm, f~14cm
    Aspher[0] = 32.0e-05 / (mm);
    Aspher[1] = -2.0e-7 / (mm3);
    Aspher[2] = 1.2e-13 / (mm2 * mm3);
  }
  else if (lens_type == 3)
  {
    Curvature = 0.0150 / mm;  // r=77mm, f~12.5cm
    Aspher[0] = 42.0e-05 / (mm);
    Aspher[1] = -3.0e-7 / (mm3);
    Aspher[2] = 1.2e-13 / (mm2 * mm3);
  }
  else if (lens_type == 4)
  {
    Curvature = 0.0175 / mm;  // r=77mm, f~10cm
    Aspher[0] = 72.0e-05 / (mm);
    Aspher[1] = -5.0e-7 / (mm3);
    Aspher[2] = 1.2e-13 / (mm2 * mm3);
  }
  else if (lens_type == 5)
  {
    Curvature = 1 / (f * (n - 1));
  }

  G4double TotAspher = 0.0 * mm;
  for (G4int k = 1; k < 9; k++)
  {
    TotAspher += Aspher[k - 1] * std::pow(r, 2 * k);
  }

  G4double ArgSqrt = 1.0 - (1.0 + Conic) * std::pow(Curvature, 2) * std::pow(r, 2);  // note conic=-1, so ArgSqrt = 1.0

  if (ArgSqrt < 0.0)
  {
    G4cout << "UltraFresnelLensParameterisation::Sagita: Square Root of <0 !" << G4endl;
  }

  G4double Sagita_value = Curvature * std::pow(r, 2) / (1.0 + std::sqrt(ArgSqrt)) + TotAspher;
  return Sagita_value;
}

//________________________________________________________________________//
PHG4mRICHDetector::mRichParameter::mRichParameter()
{
  int i;
  //----------
  // Constant
  //----------

  const double myPI = 4 * atan(1);
  fresnelLens = new LensPar();
  holderBox = new BoxPar();
  hollowVolume = new BoxPar();
  foamHolderBox = new BoxPar();
  foamHolderPoly = new PolyPar();
  aerogel = new BoxPar();
  mirror = new PolyPar();
  glassWindow = new BoxPar();
  sensor = new BoxPar();
  readout = new PolyPar();

  //--------------------Holder box key parameters-----------------------------//
  const G4double BoxDelz = 2.0 * mm;  // extract space between components
  const G4double box_thicknessXYZ[4] = {(1. / 4.) * 2.54 * cm, (1. / 4.) * 2.54 * cm, (1. / 16.) * 2.54 * cm, (1. / 4.) * 2.54 * cm};

  //------------------------- Aerogel gel key parameters---------------------//
  const G4double foamHolderThicknessXYZ[3] = {1.0 * cm, 1.0 * cm, 1.0 * cm};
  const G4double agel_halfXYZ[3] = {5.525 * cm, 5.525 * cm, 1.65 * cm};

  //------------------------Fresnel lens key parameters----------------------//
  const G4double lens_gap = (2.54 / 8.0) * cm;  //gap between agel and lens, and between lens and mirror

  const G4double lensHalfx = (5.25 / 2.0) * 2.54 * cm;
  const G4double grooveDensity = 125.0 / (2.54 * cm);

  fresnelLens->n = 1.49;
  // fresnelLens->f=6.0*2.54*cm;
  fresnelLens->f = 5.0 * 2.54 * cm;
  fresnelLens->eff_diameter = 15.24 * cm;
  fresnelLens->diameter = 2.0 * sqrt(2.0) * lensHalfx;
  fresnelLens->centerThickness = 0.06 * 2.54 * cm;
  fresnelLens->grooveWidth = (G4double) 1.0 / grooveDensity;
  fresnelLens->Set_halfXYZ(lensHalfx, grooveDensity);
  //rest of lens parameters are set below

  //---------------------------Photodetector key parameters-------------------//
  const G4double sensorGap = 0.05 * cm;  //half width of the gap
  const G4double glassWindow_halfXYZ[3] = {5.2 / 2.0 * cm, 5.2 / 2.0 * cm, 0.075 * cm};
  const G4double phodet_halfXYZ[3] = {2.4 * cm, 2.4 * cm, 0.075 * cm};

  //--------------------------- mirror key parameters ------------------------//
  const G4double mirrorThickness = 0.2 * cm;

  //-----------------------Readout Electronics key parameters-----------------//
  const G4double readout_halfz = 0.4 * cm;  //redendunt?
  const G4double readoutThickness = 0.2 * cm;

  //----------
  // calculation
  //----------
  G4double sensor_total_halfx = 2 * glassWindow_halfXYZ[0] + sensorGap;

  G4double foamHolder_halfXYZ[3];
  foamHolder_halfXYZ[0] = agel_halfXYZ[0] + foamHolderThicknessXYZ[0];
  foamHolder_halfXYZ[1] = foamHolder_halfXYZ[0];
  foamHolder_halfXYZ[2] = foamHolderThicknessXYZ[2] / 2.0;

  G4double acrylicBox_halfXYZ[3];
  acrylicBox_halfXYZ[0] = max(max(foamHolder_halfXYZ[0], sensor_total_halfx + readoutThickness), fresnelLens->halfXYZ[0]) + 0.1 * cm + box_thicknessXYZ[0];
  acrylicBox_halfXYZ[1] = acrylicBox_halfXYZ[0];
  acrylicBox_halfXYZ[2] = (BoxDelz + 2 * foamHolder_halfXYZ[2] + 2 * agel_halfXYZ[2] + lens_gap + 2 * fresnelLens->halfXYZ[2] + fresnelLens->f + 2 * glassWindow_halfXYZ[2] + 2 * phodet_halfXYZ[2] + (2 * readout_halfz + BoxDelz) + box_thicknessXYZ[2] + box_thicknessXYZ[3]) / 2.0;

  G4double hollow_halfXYZ[3];
  hollow_halfXYZ[0] = acrylicBox_halfXYZ[0] - box_thicknessXYZ[0];
  hollow_halfXYZ[1] = hollow_halfXYZ[0];
  hollow_halfXYZ[2] = (2 * acrylicBox_halfXYZ[2] - box_thicknessXYZ[2] - box_thicknessXYZ[3]) / 2.0;

  G4ThreeVector hollow_pos = G4ThreeVector(0.0 * cm, 0.0 * cm, -acrylicBox_halfXYZ[2] + hollow_halfXYZ[2] + box_thicknessXYZ[2]);

  G4double foamHolder_posz = -hollow_halfXYZ[2] + BoxDelz + foamHolder_halfXYZ[2];
  G4double agel_posz = foamHolder_posz + foamHolder_halfXYZ[2] + agel_halfXYZ[2];
  G4double lens_z = agel_posz + agel_halfXYZ[2] + fresnelLens->halfXYZ[2] + lens_gap;

  G4double glassWindow_z = lens_z - fresnelLens->halfXYZ[2] + fresnelLens->f + glassWindow_halfXYZ[2];  //out of focus. But this makes sense.
  G4double phodet_z = glassWindow_z + glassWindow_halfXYZ[2] + phodet_halfXYZ[2];

  //redendunt:
  G4double readout_z[2];
  readout_z[0] = glassWindow_z - glassWindow_halfXYZ[2];
  readout_z[1] = phodet_z + phodet_halfXYZ[2];

  //----------
  // set holderBox
  //----------
  // holderBox->name="HolderBox";
  holderBox->name = "mRICH_module";
  for (i = 0; i < 3; i++) holderBox->halfXYZ[i] = acrylicBox_halfXYZ[i];
  holderBox->pos = G4ThreeVector(0 * cm, 0 * cm, 0 * cm);
  //holderBox->material=G4Material::GetMaterial("G4_Al");
  holderBox->material = G4Material::GetMaterial("G4_AIR");  //testing
  holderBox->sensitivity = 0;

  holderBox->color = G4Colour(0.0, 0.0, 0.0);
  holderBox->visibility = true;
  holderBox->wireframe = true;
  holderBox->surface = false;

  //----------
  // set HollowVolume
  //----------
  hollowVolume->name = "HollowVolume";
  for (i = 0; i < 3; i++) hollowVolume->halfXYZ[i] = hollow_halfXYZ[i];
  hollowVolume->pos = hollow_pos;
  hollowVolume->material = G4Material::GetMaterial("mRICH_Air_Opt");
  hollowVolume->sensitivity = 0;

  hollowVolume->color = G4Colour(0.0, 0.0, 0.0);
  hollowVolume->visibility = true;
  hollowVolume->wireframe = true;
  hollowVolume->surface = false;

  //----------
  // set FoamHolder_box
  //----------
  foamHolderBox->name = "FoamHolder";
  for (i = 0; i < 3; i++) foamHolderBox->halfXYZ[i] = foamHolder_halfXYZ[i];
  foamHolderBox->pos = G4ThreeVector(0.0 * cm, 0.0 * cm, foamHolder_posz);
  foamHolderBox->material = G4Material::GetMaterial("mRICH_Air_Opt");
  foamHolderBox->sensitivity = 0;

  foamHolderBox->color = G4Colour(0.2, 0.498, 0.369);
  foamHolderBox->visibility = true;
  foamHolderBox->wireframe = true;
  foamHolderBox->surface = false;

  //----------
  // set FoamHolder_polyhedra
  //----------
  foamHolderPoly->name = "FoamHolder";
  foamHolderPoly->pos = G4ThreeVector(0, 0, 0);
  foamHolderPoly->start = 45.0 * myPI / 180.0;
  foamHolderPoly->theta = 2 * myPI;
  foamHolderPoly->numSide = 4;
  foamHolderPoly->num_zLayer = 2;

  foamHolderPoly->z[0] = agel_posz - agel_halfXYZ[2];  //front of agel
  foamHolderPoly->z[1] = agel_posz + agel_halfXYZ[2];  //back of sensor plane

  foamHolderPoly->rinner[0] = agel_halfXYZ[0];
  foamHolderPoly->rinner[1] = agel_halfXYZ[0];

  foamHolderPoly->router[0] = foamHolderPoly->rinner[0] + foamHolderThicknessXYZ[0];
  foamHolderPoly->router[1] = foamHolderPoly->rinner[1] + foamHolderThicknessXYZ[0];

  foamHolderPoly->material = G4Material::GetMaterial("mRICH_Air_Opt");
  foamHolderPoly->sensitivity = 0;

  foamHolderPoly->color = G4Colour(0.298, 0.6, 0.471);
  foamHolderPoly->visibility = true;
  foamHolderPoly->wireframe = true;
  foamHolderPoly->surface = false;

  //----------
  // set aerogel
  //----------
  aerogel->name = "Aerogel";
  for (i = 0; i < 3; i++) aerogel->halfXYZ[i] = agel_halfXYZ[i];
  aerogel->pos = G4ThreeVector(0, 0, agel_posz);
  aerogel->material = G4Material::GetMaterial("mRICH_Aerogel2");
  //aerogel->material=Air_Opt;
  aerogel->sensitivity = 0;

  aerogel->color = G4Colour(1.0, 0.65, 0.0);
  aerogel->visibility = true;
  aerogel->wireframe = true;
  aerogel->surface = false;

  //----------
  // set Fresnel lens
  //----------
  fresnelLens->name = "FresnelLens";
  fresnelLens->pos = G4ThreeVector(0, 0, lens_z);
  fresnelLens->material = G4Material::GetMaterial("mRICH_Acrylic");
  fresnelLens->sensitivity = 0;
  fresnelLens->color = G4Colour(0.0, 1.0, 1.0);
  fresnelLens->visibility = true;
  fresnelLens->wireframe = true;
  fresnelLens->surface = false;

  //----------
  // set mirror
  //----------
  mirror->name = "mirror";
  mirror->pos = G4ThreeVector(0, 0, 0);
  mirror->start = 45.0 * myPI / 180.0;
  mirror->theta = 2 * myPI;
  mirror->numSide = 4;
  mirror->num_zLayer = 2;

  mirror->z[0] = lens_z + fresnelLens->halfXYZ[2] + lens_gap;  //back of lens+air gap
  mirror->z[1] = glassWindow_z - glassWindow_halfXYZ[2];       //front of sensor plane

  mirror->rinner[0] = agel_halfXYZ[0];
  mirror->rinner[1] = sensor_total_halfx;

  mirror->router[0] = mirror->rinner[0] + mirrorThickness;
  mirror->router[1] = mirror->rinner[1] + mirrorThickness;

  mirror->material = G4Material::GetMaterial("G4_Al");
  mirror->sensitivity = 0;

  mirror->color = G4Colour(1.0, 1.0, 0.0);
  mirror->visibility = true;
  mirror->wireframe = true;

  //----------
  // set glass window
  //----------
  glassWindow->name = "glassWindow";
  for (i = 0; i < 3; i++) glassWindow->halfXYZ[i] = glassWindow_halfXYZ[i];
  glassWindow->pos = G4ThreeVector(glassWindow_halfXYZ[0] + sensorGap,  //the position of the first sensor module.
                                   glassWindow_halfXYZ[0] + sensorGap,  //the position of other sensor module will
                                   glassWindow_z);                      //be set glass window position in another func.
  glassWindow->material = G4Material::GetMaterial("mRICH_Borosilicate");
  glassWindow->sensitivity = 0;

  glassWindow->color = G4Colour(0.101, 0.737, 0.612);
  glassWindow->visibility = true;
  glassWindow->wireframe = true;
  glassWindow->surface = false;

  //----------
  // set sensor
  //----------
  sensor->name = "sensor";
  for (i = 0; i < 3; i++) sensor->halfXYZ[i] = phodet_halfXYZ[i];
  sensor->pos = G4ThreeVector(0, 0, phodet_z);  //temporary. will be set in another func.
  sensor->material = G4Material::GetMaterial("mRICH_Air_Opt");
  sensor->sensitivity = 0;

  sensor->color = G4Colour(0.0, 0.0, 0.63);
  sensor->visibility = true;
  sensor->wireframe = true;
  sensor->surface = true;

  //----------
  // set readout
  //----------
  readout->name = "readout";
  readout->pos = G4ThreeVector(0, 0, 0);
  readout->start = 45.0 * myPI / 180.0;
  readout->theta = 2 * myPI;
  readout->numSide = 4;
  readout->num_zLayer = 2;

  readout->z[0] = readout_z[0];
  readout->z[1] = readout_z[1];

  readout->rinner[0] = sensor_total_halfx;
  readout->rinner[1] = readout->rinner[0];

  readout->router[0] = readout->rinner[0] + readoutThickness;
  readout->router[1] = readout->router[0];

  readout->material = G4Material::GetMaterial("G4_Al");
  readout->sensitivity = 0;

  readout->color = G4Colour(1.0, 0.0, 0.0);
  readout->visibility = true;
  readout->wireframe = true;
  readout->surface = false;
}
//________________________________________________________________________//
PHG4mRICHDetector::mRichParameter::~mRichParameter() { ; }
//________________________________________________________________________//
void PHG4mRICHDetector::mRichParameter::SetPar_glassWindow(int i, G4double x, G4double y)
{
  glassWindow->name = "glassWindow" + std::to_string(i);
  //sprintf(glassWindow->name,"glassWindow%d",i);
  glassWindow->pos.setX(x);
  glassWindow->pos.setY(y);
}
//________________________________________________________________________//
void PHG4mRICHDetector::mRichParameter::SetPar_sensor(int i, G4double x, G4double y)
{
  sensor->name = "sensor_" + std::to_string(i);
  //sprintf(sensor->name,"sensor%d",i);
  sensor->pos.setX(x);
  sensor->pos.setY(y);
}
//________________________________________________________________________//
PHG4mRICHDetector::BoxPar* PHG4mRICHDetector::mRichParameter::GetBoxPar(string componentName)
{
  if (componentName.compare("holderBox") == 0)
    return holderBox;
  else if (componentName.compare("hollowVolume") == 0)
    return hollowVolume;
  else if (componentName.compare("foamHolderBox") == 0)
    return foamHolderBox;
  else if (componentName.compare("aerogel") == 0)
    return aerogel;
  else if (componentName.compare("glassWindow") == 0)
    return glassWindow;
  else if (componentName.compare("sensor") == 0)
    return sensor;
  else
    printf("mRichParameter::GetBoxPar() ----- ERROR: cannot find parameter=%s\n", componentName.c_str());

  return 0;
}
//________________________________________________________________________//
PHG4mRICHDetector::LensPar* PHG4mRICHDetector::mRichParameter::GetLensPar(string componentName)
{
  if (componentName.compare("fresnelLens") == 0)
    return fresnelLens;
  else
    printf("mRichParameter::GetLensPar() ----- ERROR: cannot find parameter=%s\n", componentName.c_str());
  return 0;
}
//________________________________________________________________________//
PHG4mRICHDetector::PolyPar* PHG4mRICHDetector::mRichParameter::GetPolyPar(string componentName)
{
  if (componentName.compare("foamHolderPoly") == 0)
    return foamHolderPoly;
  else if (componentName.compare("mirror") == 0)
    return mirror;
  else if (componentName.compare("readout") == 0)
    return readout;
  else
    printf("mRichParameter::GetPolyPar() ----- ERROR: cannot find parameter=%s\n", componentName.c_str());

  return 0;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_box(BoxPar* par, G4LogicalVolume* motherLV)
{
  G4Box* box = new G4Box(par->name.c_str(), par->halfXYZ[0], par->halfXYZ[1], par->halfXYZ[2]);
  G4LogicalVolume* log = new G4LogicalVolume(box, par->material, par->name.c_str(), 0, 0, 0);
  G4VPhysicalVolume* phy = new G4PVPlacement(0, par->pos, log, par->name.c_str(), motherLV, false, 0, OverlapCheck());

  G4VisAttributes* visAtt = new G4VisAttributes(par->color);
  visAtt->SetVisibility(par->visibility);
  visAtt->SetForceWireframe(par->wireframe);
  visAtt->SetForceSolid(par->surface);
  log->SetVisAttributes(visAtt);

  return phy;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_polyhedra(PolyPar* par, G4LogicalVolume* motherLV)
{
  G4Polyhedra* polyhedra = new G4Polyhedra(par->name.c_str(), par->start, par->theta, par->numSide,
                                           par->num_zLayer, par->z, par->rinner, par->router);
  G4LogicalVolume* log = new G4LogicalVolume(polyhedra, par->material, par->name.c_str(), 0, 0, 0);
  G4VPhysicalVolume* phy = new G4PVPlacement(0, par->pos, log, par->name.c_str(), motherLV, false, 0, OverlapCheck());

  G4VisAttributes* visAtt = new G4VisAttributes(par->color);
  visAtt->SetVisibility(par->visibility);
  visAtt->SetForceWireframe(par->wireframe);
  visAtt->SetForceSolid(par->surface);
  log->SetVisAttributes(visAtt);

  return phy;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_holderBox(mRichParameter* detectorParameter, G4LogicalVolume* motherLV)
{
  mRICH_PV = build_box(detectorParameter->GetBoxPar("holderBox"), motherLV);
  return build_box(detectorParameter->GetBoxPar("hollowVolume"), mRICH_PV->GetLogicalVolume());
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_foamHolder(mRichParameter* detectorParameter, G4LogicalVolume* motherLV)
{
  G4VPhysicalVolume* box = build_box(detectorParameter->GetBoxPar("foamHolderBox"), motherLV);
  if (!box)
  {
    cout << "placement of foamholderbox failed" << endl;
  }
  box = build_polyhedra(detectorParameter->GetPolyPar("foamHolderPoly"), motherLV);
  if (!box)
  {
    cout << "placement of foamholderpoly failed" << endl;
  }
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_aerogel(mRichParameter* detectorParameter, G4VPhysicalVolume* motherPV)
{
  G4VPhysicalVolume* aerogel = build_box(detectorParameter->GetBoxPar("aerogel"), motherPV->GetLogicalVolume());
  aerogel_vol[aerogel] = 0;

  const G4int num = 2;
  G4double Ephoton[num] = {2.034 * eV, 4.136 * eV};

  G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
  OpWaterSurface->SetType(dielectric_dielectric);
  OpWaterSurface->SetFinish(ground);
  OpWaterSurface->SetModel(unified);
  new G4LogicalBorderSurface("WaterSurface", aerogel, motherPV, OpWaterSurface);

  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num] = {0.3, 0.3};
  G4double SpecularSpike[num] = {0.2, 0.2};
  G4double Backscatter[num] = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  myST1->AddProperty("RINDEX", Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT", Ephoton, SpecularLobe, num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike, num);
  myST1->AddProperty("BACKSCATTERCONSTANT", Ephoton, Backscatter, num);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_mirror(mRichParameter* detectorParameter, G4VPhysicalVolume* motherPV)
{
  G4VPhysicalVolume* mirror = build_polyhedra(detectorParameter->GetPolyPar("mirror"), motherPV->GetLogicalVolume());

  //-----------
  //   Optical properties of the interface between the Air and Reflective Surface
  //   For Mirror, reflectivity is set at 95% and specular reflection is assumed.
  //-----------
  G4OpticalSurface* OpticalAirMirror = new G4OpticalSurface("AirMirrorSurface");
  OpticalAirMirror->SetModel(unified);
  OpticalAirMirror->SetType(dielectric_dielectric);
  OpticalAirMirror->SetFinish(polishedfrontpainted);

  const G4int NUM = 2;
  G4double lambda_min = 200 * nm;
  G4double lambda_max = 700 * nm;

  G4double XX[NUM] = {h_Planck * c_light / lambda_max, h_Planck * c_light / lambda_min};
  G4double ICEREFLECTIVITY[NUM] = {0.95, 0.95};

  G4MaterialPropertiesTable* AirMirrorMPT = new G4MaterialPropertiesTable();
  AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY, NUM);
  OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);

  new G4LogicalBorderSurface("Air/Mirror Surface", motherPV, mirror, OpticalAirMirror);
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_sensor(mRichParameter* detectorParameter, G4LogicalVolume* motherLV)
{
  //position of the first sensor module
  G4double last_x = detectorParameter->GetBoxPar("glassWindow")->pos.getX();
  G4double last_y = last_x;

  G4double x, y;

  int i;
  for (i = 0; i < 4; i++)
  {
    if (i == 0)
    {
      x = last_x;
      y = last_y;
    }
    else
    {
      x = -last_y;
      y = last_x;
    }

    detectorParameter->SetPar_glassWindow(i + 1, x, y);

    detectorParameter->SetPar_sensor(i + 1, x, y);
    sensor_PV[i] = build_box(detectorParameter->GetBoxPar("sensor"), motherLV);

    sensor_vol[sensor_PV[i]] = i;
    // cout << "in build_sensor: sensor_vol = " << sensor_vol[sensor_PV[i]] << endl;

    last_x = x;
    last_y = y;
  }  //end of for(i)
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_lens(LensPar* par, G4LogicalVolume* motherLV)
{
  const G4int NumberOfGrooves = floor((par->eff_diameter / 2.0) / par->grooveWidth);
  G4Polycone* Groove_poly[NumberOfGrooves];
  G4LogicalVolume* Groove_log[NumberOfGrooves];

  G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  SurfaceVisAtt->SetVisibility(true);
  SurfaceVisAtt->SetForceWireframe(true);
  SurfaceVisAtt->SetForceSolid(true);

  int igroove;
  for (igroove = 0; igroove < 1000; igroove++)
  {  //just put a arbitrary large number
    //--------------------------------
    //Grooves' inner and outer radius
    //--------------------------------
    G4double iRmin1 = (igroove + 0) * par->grooveWidth;
    G4double iRmax1 = (igroove + 1) * par->grooveWidth;
    G4double iRmin2 = iRmin1;
    G4double iRmax2 = iRmin2 + 0.0001;

    G4double lens_poly_rmin[3] = {iRmin1, iRmin1, iRmin2};
    G4double lens_poly_rmax[3] = {iRmax1, iRmax1, iRmax2};

    if (iRmax1 > par->diameter / 2.0) break;  //if iRmax1>Lens radius (outside the lens), break

    //--------------------------------
    //phi angle
    //--------------------------------
    G4double phi1;
    G4double phi2;
    G4double deltaPhi;

    if (iRmax1 < par->halfXYZ[0])
    {                    //draw a full circle
      phi1 = 0;          //in rad
      deltaPhi = twopi;  //in rad. two pi
    }
    else
    {
      phi1 = acos(par->halfXYZ[0] / iRmax1);  //in rad
      phi2 = asin(par->halfXYZ[0] / iRmax1);  //in rad, assume lens is square -> halfy=halfx
      deltaPhi = phi2 - phi1;
    }
    //--------------------------------
    //grooves profile
    //--------------------------------
    G4double lens_poly_z[3];
    int numOfLayer;

    if (iRmin1 < par->eff_diameter / 2.0)
    {  //if iRmin>=effective radius, dZ=0, i.e. flat
      numOfLayer = 3;
      G4double dZ = par->GetSagita(iRmax1) - par->GetSagita(iRmin1);
      lens_poly_z[0] = par->halfXYZ[2];
      lens_poly_z[1] = -par->halfXYZ[2] + dZ;
      lens_poly_z[2] = -par->halfXYZ[2];
    }
    else
    {
      numOfLayer = 2;
      lens_poly_z[0] = par->halfXYZ[2];
      lens_poly_z[1] = par->halfXYZ[2] - par->centerThickness;
      lens_poly_z[2] = 0;
    }

    //--------------------------------
    //build grooves
    //--------------------------------
    int repeat = 1;
    if (iRmax1 >= par->halfXYZ[0])
    {
      repeat = 4;
    }  //4 edges
    for (int i = 0; i < repeat; i++)
    {
      Groove_poly[i] = new G4Polycone(par->name.c_str(), phi1, deltaPhi, numOfLayer, lens_poly_z, lens_poly_rmin, lens_poly_rmax);
      Groove_log[i] = new G4LogicalVolume(Groove_poly[i], par->material, par->name.c_str(), 0, 0, 0);
      new G4PVPlacement(0, par->pos, Groove_log[i], par->name.c_str(), motherLV, false, 0, OverlapCheck());

      Groove_log[i]->SetVisAttributes(SurfaceVisAtt);
      phi1 = phi1 + halfpi;  //g4 pre-defined: halfpi=pi/2
    }
  }
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_mRICH_wall_hside(G4LogicalVolume* logicWorld)
{
  G4AssemblyVolume* mRICHwall = new G4AssemblyVolume();  //"mother volume"

  G4LogicalVolume* a_mRICH = Construct_a_mRICH(0);  // build a single mRICH

  int NumOfModule = params->get_int_param("NumOfModule_wall_hside");

  for (int i_mRICH = 0; i_mRICH < NumOfModule; ++i_mRICH)
  {
    // get moduleID
    // std::stringstream key_moduleID;
    // key_moduleID << "mRICH_wall_hside_" << i_mRICH << "_moduleID";
    // int module_id = params->get_int_param(key_moduleID.str());

    // get position
    std::stringstream key_position_x;
    key_position_x << "mRICH_wall_hside_" << i_mRICH << "_position_x";
    G4double x = params->get_double_param(key_position_x.str());

    std::stringstream key_position_y;
    key_position_y << "mRICH_wall_hside_" << i_mRICH << "_position_y";
    G4double y = params->get_double_param(key_position_y.str());

    std::stringstream key_position_z;
    key_position_z << "mRICH_wall_hside_" << i_mRICH << "_position_z";
    G4double z = params->get_double_param(key_position_z.str());

    // get rotation
    std::stringstream key_rotation_theta;
    key_rotation_theta << "mRICH_wall_hside_" << i_mRICH << "_rotation_theta";
    G4double theta = params->get_double_param(key_rotation_theta.str());

    std::stringstream key_rotation_phi;
    key_rotation_phi << "mRICH_wall_hside_" << i_mRICH << "_rotation_phi";
    G4double phi = params->get_double_param(key_rotation_phi.str());

    // cout << "module_id = " << module_id << ", x = " << x << ", y = " << y << ", z = " << z << ", theta = " << theta << ", phi = " << phi << endl;

    G4ThreeVector pos(x, y, z);
    G4RotationMatrix* rot = new G4RotationMatrix();

    if (x != 0 || y != 0)
    {
      rot->rotateX(theta * (-1) * sin(phi) * 180 * deg / pi);
      rot->rotateY(theta * cos(phi) * 180 * deg / pi);
    }
    mRICHwall->AddPlacedVolume(a_mRICH, pos, rot);
  }

  G4ThreeVector pos(0, 0, 0);
  mRICHwall->MakeImprint(logicWorld, pos, nullptr, 0, OverlapCheck());

  printf("-----------------------------------------------------------------------------\n");
  printf("%d detectors are built\n", NumOfModule);
  printf("-----------------------------------------------------------------------------\n");
}
//________________________________________________________________________//

//________________________________________________________________________//
void PHG4mRICHDetector::build_mRICH_wall_eside(G4LogicalVolume* logicWorld)
{
  G4AssemblyVolume* mRICHwall = new G4AssemblyVolume();  //"mother volume"

  G4LogicalVolume* a_mRICH = Construct_a_mRICH(0);  // build a single mRICH

  G4double shift = params->get_double_param("mRICH_wall_eside_shift");

  int NumOfModule = params->get_int_param("NumOfModule_wall_eside");

  for (int i_mRICH = 0; i_mRICH < NumOfModule; ++i_mRICH)
  {
    // get moduleID
    // std::stringstream key_moduleID;
    // key_moduleID << "mRICH_wall_eside_" << i_mRICH << "_moduleID";
    // int module_id = params->get_int_param(key_moduleID.str());

    // get position
    std::stringstream key_position_x;
    key_position_x << "mRICH_wall_eside_" << i_mRICH << "_position_x";
    G4double x = params->get_double_param(key_position_x.str());

    std::stringstream key_position_y;
    key_position_y << "mRICH_wall_eside_" << i_mRICH << "_position_y";
    G4double y = params->get_double_param(key_position_y.str());

    std::stringstream key_position_z;
    key_position_z << "mRICH_wall_eside_" << i_mRICH << "_position_z";
    G4double z = params->get_double_param(key_position_z.str());

    // cout << "module_id = " << module_id << ", x = " << x << ", y = " << y << ", z = " << z << endl;

    G4ThreeVector pos(x, y, z);
    G4RotationMatrix* rot = new G4RotationMatrix();
    mRICHwall->AddPlacedVolume(a_mRICH, pos, rot);
  }

  G4ThreeVector pos(0, 0, shift);
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot->rotateX(180 * deg);
  mRICHwall->MakeImprint(logicWorld, pos, rot, 0, OverlapCheck());

  printf("-----------------------------------------------------------------------------\n");
  printf("%d detectors are built\n", NumOfModule);
  printf("-----------------------------------------------------------------------------\n");
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_mRICH_sector(G4LogicalVolume* logicWorld, int numSector)
{
  G4AssemblyVolume* sector = new G4AssemblyVolume();  //"mother volume"

  G4LogicalVolume* a_mRICH = Construct_a_mRICH(0);  // build a single mRICH

  G4double theta = params->get_double_param("mRICH_sector_hside_rotation_theta");

  G4double shift = params->get_double_param("mRICH_sector_hside_shift");

  int NumOfModule = params->get_int_param("NumOfModule_sector_hside");

  for (int i_mRICH = 0; i_mRICH < NumOfModule; ++i_mRICH)
  {
    // get moduleID
    // std::stringstream key_moduleID;
    // key_moduleID << "mRICH_sector_hside_" << i_mRICH << "_moduleID";
    // int module_id = params->get_int_param(key_moduleID.str());

    // get position
    std::stringstream key_position_x;
    key_position_x << "mRICH_sector_hside_" << i_mRICH << "_position_x";
    G4double x = params->get_double_param(key_position_x.str());

    std::stringstream key_position_y;
    key_position_y << "mRICH_sector_hside_" << i_mRICH << "_position_y";
    G4double y = params->get_double_param(key_position_y.str());

    std::stringstream key_position_z;
    key_position_z << "mRICH_sector_hside_" << i_mRICH << "_position_z";
    G4double z = params->get_double_param(key_position_z.str());

    // cout << "module_id = " << module_id << ", x = " << x << ", y = " << y << ", z = " << z << ", theta = " << theta << endl;

    G4ThreeVector pos(x, y, z);
    G4RotationMatrix* rot = new G4RotationMatrix();

    sector->AddPlacedVolume(a_mRICH, pos, rot);
  }

  for (int i = 0; i < numSector; i++)
  {
    // G4ThreeVector pos(0, 0, 3.0*m);
    G4ThreeVector pos(0, 0, shift);
    // G4ThreeVector pos(0, 0, 2.8085*m);
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateX(-theta * 180 * deg / pi);
    rot->rotateZ(i * 45 * deg);
    sector->MakeImprint(logicWorld, pos, rot, 0, OverlapCheck());
  }
}

//________________________________________________________________________//
