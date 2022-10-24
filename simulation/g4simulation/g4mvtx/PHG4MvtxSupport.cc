#include "PHG4MvtxSupport.h"

#include "PHG4MvtxCable.h"
#include "PHG4MvtxDefs.h"
#include "PHG4MvtxDisplayAction.h"
#include "PHG4MvtxServiceStructure.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4UserLimits.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/format.hpp>
#pragma GCC diagnostic pop

#include <algorithm>  // for max
#include <cmath>      // for M_PI, atan, cos, sin
#include <cstddef>    // for NULL
#include <iostream>   // for operator<<, basic_...
#include <utility>    // for pair, make_pair

class G4VSolid;

namespace ServiceProperties
{
  std::string materials[] = {"G4_Cu", "MVTX_CarbonFiber$", "G4_POLYETHYLENE"};
  const int nMaterials = sizeof(materials) / sizeof(materials[0]);

  float ServiceEnd = -7.21;
  float ServiceOffset = -15.0;
  float BarrelOffset = 18.679;
  float BarrelRadius = 10.33;     //Inner radious of service barrel
  float BarrelThickness = 0.436;  //Thickness in cm
  float BarrelLength = 121.24;    //Length of cylinder in cm
  float BarrelCableStart = -1. * BarrelOffset - 25.;
  float LayerThickness = 0.1;  //
  float CYSSConeThickness = 0.216;
  float CYSSRibThickness = 0.170;
  float cableRotate[3] = {10., 5., 5.};  //Rotate the cables to line up with the staves
  float radToDeg = 180.0 / M_PI;
  float degToRad = 1. / radToDeg;
}  // namespace ServiceProperties

using namespace ServiceProperties;

PHG4MvtxSupport::PHG4MvtxSupport(PHG4MvtxDisplayAction *dispAct, bool overlapCheck)
  : m_DisplayAction(dispAct)
  , m_avSupport(nullptr)
  , m_avBarrelCable(nullptr)
  , m_avL0Cable(nullptr)
  , m_avL1Cable(nullptr)
  , m_avL2Cable(nullptr)
{
  m_overlapCheck = overlapCheck;
}

PHG4MvtxSupport::~PHG4MvtxSupport()
{
  delete m_avSupport;
  delete m_avBarrelCable;
  delete m_avL0Cable;
  delete m_avL1Cable;
  delete m_avL2Cable;
}

std::vector<float> PHG4MvtxSupport::get_thickness(PHG4MvtxServiceStructure *object)
{
  std::vector<float> thickness = {object->get_thickness_copper(), object->get_thickness_carbon(), object->get_thickness_plastic()};
  return thickness;
}

G4Material *supportMaterial()
{
  G4double density;
  G4int natoms;
  G4String symbol;

  G4Material *Epoxy = new G4Material("MVTX_Epoxy$", density = 1.56 * g / cm3, natoms = 4);
  Epoxy->AddElement(PHG4Detector::GetDetectorElement("H", true), 32);  // Hydrogen
  Epoxy->AddElement(PHG4Detector::GetDetectorElement("C", true), 2);   // Nitrogen
  Epoxy->AddElement(PHG4Detector::GetDetectorElement("N", true), 4);   // Oxygen
  Epoxy->AddElement(PHG4Detector::GetDetectorElement("O", true), 15);  // Carbon

  G4Material *G4_mat = new G4Material("MVTX_CarbonFiber$", density = 1.987 * g / cm3, natoms = 2);
  G4_mat->AddMaterial(PHG4Detector::GetDetectorMaterial("G4_C"), 70 * perCent);  // Carbon (NX-80-240)
  G4_mat->AddMaterial(Epoxy, 30 * perCent);                                      // Epoxy (EX-1515)

  return G4_mat;
}

G4Material *carbonFibreMaterial = supportMaterial();

void PHG4MvtxSupport::TrackingServiceCone(PHG4MvtxServiceStructure *object, G4AssemblyVolume &assemblyVolume)
{
  float length = std::abs(object->get_zNorth() - object->get_zSouth());
  std::vector<float> thickness = get_thickness(object);
  float innerRadiusSouth = object->get_rSouth();
  float innerRadiusNorth = object->get_rNorth();
  float outerRadiusSouth;
  float outerRadiusNorth;

  G4RotationMatrix *rot = new G4RotationMatrix();
  rot->rotateZ(0.);
  G4ThreeVector place;
  place.setZ((object->get_zSouth() + length / 2) * cm);

  for (int i = 0; i < nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    outerRadiusSouth = innerRadiusSouth + thickness[i];
    outerRadiusNorth = innerRadiusNorth + thickness[i];

    G4Material *trackerMaterial = nullptr;
    if (materials[i] == "MVTX_CarbonFiber$")
      trackerMaterial = carbonFibreMaterial;
    else
      trackerMaterial = PHG4Detector::GetDetectorMaterial(materials[i]);

    G4VSolid *coneSolid = new G4Cons(G4String(object->get_name() + "_SOLID"),
                                     innerRadiusSouth * cm, outerRadiusSouth * cm,
                                     innerRadiusNorth * cm, outerRadiusNorth * cm, (length / 2) * cm, 0, 2 * M_PI);

    G4LogicalVolume *coneLogic = new G4LogicalVolume(coneSolid, trackerMaterial,
                                                     G4String(object->get_name() + "_LOGIC"), nullptr, nullptr, nullptr);

    m_DisplayAction->AddVolume(coneLogic, materials[i]);

    assemblyVolume.AddPlacedVolume(coneLogic, place, rot);

    innerRadiusSouth = outerRadiusSouth;
    innerRadiusNorth = outerRadiusNorth;
  }
}

void PHG4MvtxSupport::TrackingServiceCylinder(PHG4MvtxServiceStructure *object, G4AssemblyVolume &assemblyVolume)
{
  float length = std::abs(object->get_zNorth() - object->get_zSouth());
  std::vector<float> thickness = get_thickness(object);
  float innerRadius = object->get_rSouth();
  float outerRadius;

  G4RotationMatrix *rot = new G4RotationMatrix();
  rot->rotateZ(0.);
  G4ThreeVector place;
  place.setZ((object->get_zSouth() + length / 2) * cm);

  for (int i = 0; i < nMaterials; ++i)
  {
    if (thickness[i] == 0) continue;
    outerRadius = innerRadius + thickness[i];

    G4Material *trackerMaterial = nullptr;
    if (materials[i] == "MVTX_CarbonFiber$")
    {
      trackerMaterial = carbonFibreMaterial;
    }
    else
    {
      trackerMaterial = PHG4Detector::GetDetectorMaterial(materials[i]);
    }
    G4VSolid *cylinderSolid = new G4Tubs(G4String(object->get_name() + "_SOLID"),
                                         innerRadius * cm, outerRadius * cm, (length / 2) * cm, 0, 2 * M_PI);

    G4LogicalVolume *cylinderLogic = new G4LogicalVolume(cylinderSolid, trackerMaterial,
                                                         G4String(object->get_name() + "_LOGIC"), nullptr, nullptr, nullptr);

    m_DisplayAction->AddVolume(cylinderLogic, materials[i]);

    assemblyVolume.AddPlacedVolume(cylinderLogic, place, rot);

    innerRadius = outerRadius;
  }
}

void PHG4MvtxSupport::CreateCable(PHG4MvtxCable *object, G4AssemblyVolume &assemblyVolume)
{
  std::string cableMaterials[2] = {object->get_coreMaterial(), "G4_POLYETHYLENE"};

  float dX = object->get_xNorth() - object->get_xSouth();
  float dY = object->get_yNorth() - object->get_ySouth();
  float dZ = object->get_zNorth() - object->get_zSouth();

  float rotY = dZ != 0. ? std::atan(dX / dZ) : 0.;
  float rotZ = dX != 0. ? std::atan(dY / dX) : 0.;

  float setX = (object->get_xSouth() + object->get_xNorth()) / 2;
  float setY = (object->get_ySouth() + object->get_yNorth()) / 2;
  float setZ = (object->get_zSouth() + object->get_zNorth()) / 2;

  float length = std::sqrt(dX * dX + dY * dY + dZ * dZ);
  float IR[2] = {0, object->get_coreRadius()};
  float OR[2] = {object->get_coreRadius(), object->get_sheathRadius()};

  G4RotationMatrix rot;
  rot.rotateY(rotY);
  rot.rotateZ(rotZ);
  G4ThreeVector place;
  place.setX(setX * cm);
  place.setY(setY * cm);
  place.setZ(setZ * cm);
  G4Transform3D transform(rot, place);
  // we need just one of these but have multiple calls to this method
  static G4UserLimits *g4userLimits = new G4UserLimits(0.01);

  for (int i = 0; i < 2; ++i)
  {
    G4Material *trackerMaterial = PHG4Detector::GetDetectorMaterial(cableMaterials[i]);

    G4VSolid *cylinderSolid = new G4Tubs(G4String(object->get_name() + "_SOLID"),
                                         IR[i] * cm, OR[i] * cm, (length / 2.) * cm, 0, 2 * M_PI);

    G4LogicalVolume *cylinderLogic = new G4LogicalVolume(cylinderSolid, trackerMaterial,
                                                         G4String(object->get_name() + "_LOGIC"), nullptr, nullptr, g4userLimits);

    if (i == 0)
    {
      m_DisplayAction->AddVolume(cylinderLogic, cableMaterials[i]);
    }
    else
    {
      m_DisplayAction->AddVolume(cylinderLogic, object->get_color());
    }

    assemblyVolume.AddPlacedVolume(cylinderLogic, transform);
  }
}

void PHG4MvtxSupport::CreateCableBundle(G4AssemblyVolume &assemblyVolume, const std::string &superName,
                                        bool enableSignal, bool enableCooling, bool enablePower,
                                        float x1, float x2, float y1, float y2, float z1, float z2)  //, float theta)
{
  //Set up basic MVTX cable bundle (24 Samtec cables, 1 power cable, 2 cooling cables)
  float samtecCoreRadius = 0.01275;
  float samtecSheathRadius = 0.05;
  float coolingCoreRadius = 0.056;
  float coolingSheathRadius = 0.2;  //?
  float powerLargeCoreRadius = 0.069;
  float powerLargeSheathRadius = 0.158;
  float powerMediumCoreRadius = 0.033;
  float powerMediumSheathRadius = 0.082;
  float powerSmallCoreRadius = 0.028;
  float powerSmallSheathRadius = 0.0573;  //?

  float globalShiftX = 0.;
  float globalShiftY = -0.0984;
  float samtecShiftX = -6 * samtecSheathRadius + globalShiftX;
  float samtecShiftY = 1 * samtecSheathRadius + globalShiftY;
  float coolingShiftX = -3 * coolingSheathRadius + globalShiftX;
  float coolingShiftY = -1 * coolingSheathRadius + globalShiftY;
  float powerShiftX = 3.5 * powerLargeSheathRadius + globalShiftX;
  float powerShiftY = 6.1 * powerLargeSheathRadius + globalShiftY;

  //Samtec cables (we use 24 as there are 12 twinax)
  if (enableSignal)
  {
    unsigned int nSamtecWires = 24;
    unsigned int nRow = 8;
    unsigned int nCol = nSamtecWires / nRow;
    for (unsigned int iRow = 0; iRow < nRow; ++iRow)
    {
      for (unsigned int iCol = 0; iCol < nCol; ++iCol)
      {
        float deltaX = samtecShiftX + ((iCol + 1) * (samtecSheathRadius * 2.6));
        float deltaY = samtecShiftY - ((iRow + 1) * (samtecSheathRadius * 2.1));
        PHG4MvtxCable *cable = new PHG4MvtxCable(boost::str(boost::format("%s_samtec_%d_%d") % superName.c_str() % iRow % iCol), "G4_Cu", samtecCoreRadius, samtecSheathRadius,
                                                 x1 + deltaX, x2 + deltaX, y1 + deltaY, y2 + deltaY, z1, z2, "blue");
        CreateCable(cable, assemblyVolume);
        delete cable;
      }
    }
  }

  //Cooling Cables
  if (enableCooling)
  {
    unsigned int nCool = 2;
    std::string cooling_color[2] = {"red", "white"};
    for (unsigned int iCool = 0; iCool < nCool; ++iCool)
    {
      float deltaX = coolingShiftX + ((iCool + 1) * (coolingSheathRadius * 2));
      float deltaY = coolingShiftY + (coolingSheathRadius * 2);
      PHG4MvtxCable *cable = new PHG4MvtxCable(boost::str(boost::format("%s_cooling_%d") % superName.c_str() % iCool), "G4_WATER", coolingCoreRadius, coolingSheathRadius,
                                               x1 + deltaX, x2 + deltaX, y1 + deltaY, y2 + deltaY, z1, z2, cooling_color[iCool]);
      CreateCable(cable, assemblyVolume);
      delete cable;
    }
  }

  //Power Cables
  if (enablePower)
  {
    using PowerCableParameters = std::pair<std::pair<std::string, std::string>, std::pair<float, float>>;
    std::vector<PowerCableParameters> powerCables;
    std::vector<std::vector<float>> powerCableColors;

    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiReturn") % superName.c_str()), "Large"), std::make_pair((-2.5 * powerLargeSheathRadius) + powerShiftX, (-2.5 * powerLargeSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiSupply") % superName.c_str()), "Large"), std::make_pair((-4.5 * powerLargeSheathRadius) + powerShiftX, (-1.5 * powerLargeSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaReturn") % superName.c_str()), "Medium"), std::make_pair((-4 * powerLargeSheathRadius) + powerShiftX, (-5.75 * powerMediumSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaSupply") % superName.c_str()), "Medium"), std::make_pair((-5.1 * powerLargeSheathRadius) + powerShiftX, (-5.75 * powerMediumSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_digiSense") % superName.c_str()), "Small"), std::make_pair((-10 * powerSmallSheathRadius) + powerShiftX, (-1 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_anaSense") % superName.c_str()), "Small"), std::make_pair((-8 * powerSmallSheathRadius) + powerShiftX, (-2 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_bias") % superName.c_str()), "Small"), std::make_pair((-6 * powerSmallSheathRadius) + powerShiftX, (-3 * powerSmallSheathRadius) + powerShiftY));
    powerCables.emplace_back(std::make_pair(boost::str(boost::format("%s_ground") % superName.c_str()), "Small"), std::make_pair((-4 * powerSmallSheathRadius) + powerShiftX, (-4 * powerSmallSheathRadius) + powerShiftY));

    for (PowerCableParameters &powerCable : powerCables)
    {
      float coreRad, sheathRad;
      std::string cableColor;
      std::string cableType = powerCable.first.second;
      std::string cableName = powerCable.first.first;
      if (cableType == "Small")
      {
        coreRad = powerSmallCoreRadius;
        sheathRad = powerSmallSheathRadius;
      }
      else if (cableType == "Medium")
      {
        coreRad = powerMediumCoreRadius;
        sheathRad = powerMediumSheathRadius;
      }
      else
      {
        coreRad = powerLargeCoreRadius;
        sheathRad = powerLargeSheathRadius;
      }

      if (cableName == boost::str(boost::format("%s_digiReturn") % superName.c_str())) cableColor = "black";
      if (cableName == boost::str(boost::format("%s_digiSupply") % superName.c_str())) cableColor = "red";
      if (cableName == boost::str(boost::format("%s_anaReturn") % superName.c_str())) cableColor = "black";
      if (cableName == boost::str(boost::format("%s_anaSupply") % superName.c_str())) cableColor = "red";
      if (cableName == boost::str(boost::format("%s_digiSense") % superName.c_str())) cableColor = "white";
      if (cableName == boost::str(boost::format("%s_anaSense") % superName.c_str())) cableColor = "green";
      if (cableName == boost::str(boost::format("%s_bias") % superName.c_str())) cableColor = "white";
      if (cableName == boost::str(boost::format("%s_ground") % superName.c_str())) cableColor = "green";

      PHG4MvtxCable *cable = new PHG4MvtxCable(powerCable.first.first, "G4_Cu", coreRad, sheathRad,
                                               (x1 + powerCable.second.first), (x2 + powerCable.second.first),
                                               (y1 + powerCable.second.second), (y2 + powerCable.second.second), z1, z2, cableColor);
      CreateCable(cable, assemblyVolume);
      delete cable;
    }
  }
}

G4AssemblyVolume *PHG4MvtxSupport::buildBarrelCable()
{
  G4AssemblyVolume *av = new G4AssemblyVolume();

  CreateCableBundle(*av, "barrelCable", true, true, true, 0, 0, 0, 0, -1. * (BarrelLength + BarrelOffset), BarrelCableStart);

  return av;
}

G4AssemblyVolume *PHG4MvtxSupport::buildL0Cable()
{
  G4AssemblyVolume *av = new G4AssemblyVolume();
  float rInner = 2.297;
  float rOuter = 4.250;
  float zMin = -18.680;
  float zTransition1 = -17.079;
  float zTransition2 = -9.186;
  float zMax = ServiceEnd;
  CreateCableBundle(*av, "MVTX_L0Cable_0", true, false, false, rOuter, rOuter, 0, 0, zMin, zTransition1);
  CreateCableBundle(*av, "MVTX_L0Cable_1", true, false, false, rOuter, rInner, 0, 0, zTransition1 + 0.1, zTransition2);
  CreateCableBundle(*av, "MVTX_L0Cable_2", true, false, false, rInner, rInner, 0, 0, zTransition2 + 0.1, zMax);

  return av;
}

G4AssemblyVolume *PHG4MvtxSupport::buildL1Cable()
{
  G4AssemblyVolume *av = new G4AssemblyVolume();
  float rInner = 3.299;
  float rOuter = 6.738;
  float zMin = -18.000;
  float zTransition1 = -15.851;
  float zTransition2 = -8.938;
  float zMax = ServiceEnd;
  CreateCableBundle(*av, "MVTX_L1Cable_0", true, false, false, rOuter, rOuter, 0, 0, zMin, zTransition1);
  CreateCableBundle(*av, "MVTX_L1Cable_1", true, false, false, rOuter, rInner, 0, 0, zTransition1 + 0.1, zTransition2);
  CreateCableBundle(*av, "MVTX_L1Cable_2", true, false, false, rInner, rInner, 0, 0, zTransition2 + 0.1, zMax);

  return av;
}

G4AssemblyVolume *PHG4MvtxSupport::buildL2Cable()
{
  G4AssemblyVolume *av = new G4AssemblyVolume();
  float rInner = 4.074;
  float rOuter = 9.080;
  float zMin = -22.300;
  float zTransition1 = -15.206;
  float zTransition2 = -8.538;
  float zMax = ServiceEnd;
  CreateCableBundle(*av, "MVTX_L2Cable_0", true, false, false, rOuter, rOuter, 0, 0, zMin, zTransition1);
  CreateCableBundle(*av, "MVTX_L2Cable_1", true, false, false, rOuter, rInner, 0, 0, zTransition1 + 0.1, zTransition2);
  CreateCableBundle(*av, "MVTX_L2Cable_2", true, false, false, rInner, rInner, 0, 0, zTransition2 + 0.1, zMax);

  return av;
}

void PHG4MvtxSupport::ConstructMvtxSupport(G4LogicalVolume *&lv)
{
  unsigned int nStaves[PHG4MvtxDefs::kNLayers];
  unsigned int totStaves = 0;
  for (unsigned int i = 0; i < PHG4MvtxDefs::kNLayers; ++i)
  {
    nStaves[i] = (int) PHG4MvtxDefs::mvtxdat[i][PHG4MvtxDefs::kNStave];
    totStaves += nStaves[i];
  }

  std::vector<PHG4MvtxServiceStructure *> cylinders, cones;
  m_avSupport = new G4AssemblyVolume();

  //Service Barrel
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_serviceBarrel_0", 0, BarrelThickness, 0, -1. * (BarrelLength + BarrelOffset), -26.209, BarrelRadius, 0));

  //CYSS
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Cone_0", 0, CYSSConeThickness, 0., -26.208, -15.68, 10.55, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Cone_1", 0, CYSSConeThickness, 0., -15.679, -8.619, 10.55, 5.302));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Cone_2", 0, CYSSConeThickness, 0., -8.618, -6.18, 5.302, 0));

  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Rib_0", 0, CYSSRibThickness, 0., -21.719, -20.949, 9.762, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Rib_1", 0, CYSSRibThickness, 0., -20.948, -20.159, 9.762, 10.36));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Rib_2", 0, CYSSRibThickness, 0., -20.158, -17.749, 10.36, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Rib_3", 0, CYSSRibThickness, 0., -17.748, -16.959, 10.36, 9.762));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Rib_4", 0, CYSSRibThickness, 0., -16.958, -16.196, 9.762, 0));

  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_CYSS_Cylinder", 0, 0.112, 0, -8.619, 36.153, 5.15, 0));

  //MVTX Layers
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L0_0", 0, LayerThickness, 0., -18.680, -16.579, 5.050, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_L0_1", 0, LayerThickness, 0., -16.578, -9.186, 5.050, 2.997));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L0_2", 0, LayerThickness, 0., -9.185, ServiceEnd, 2.997, 0));

  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L1_0", 0, LayerThickness, 0., -17.970, -15.851, 7.338, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_L1_1", 0, LayerThickness, 0., -15.850, -8.938, 7.338, 3.799));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L1_2", 0, LayerThickness, 0., -8.937, ServiceEnd, 3.799, 0));

  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L2_0", 0, LayerThickness, 0., -22.300, -15.206, 9.580, 0));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_L2_1", 0, LayerThickness, 0., -15.205, -8.538, 9.580, 4.574));
  cylinders.push_back(new PHG4MvtxServiceStructure("MVTX_L2_2", 0, LayerThickness, 0., -8.537, ServiceEnd, 4.574, 0));

  //Conenct copper from barrel to layers
  //Currently non-discrete cones as rotations are acting up
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_sb_to_L0", 0.005, 0., 0.066, -26.9, -18.680, 10.10, 5.050));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_sb_to_L1", 0.004, 0., 0.061, -26.9, -18.000, 10.20, 7.338));
  cones.push_back(new PHG4MvtxServiceStructure("MVTX_sb_to_L2", 0.004, 0., 0.058, -26.7, -22.301, 10.25, 9.580));

  for (PHG4MvtxServiceStructure *cylinder : cylinders)
  {
    TrackingServiceCylinder(cylinder, *m_avSupport);
    delete cylinder;
  }
  for (PHG4MvtxServiceStructure *cone : cones)
  {
    TrackingServiceCone(cone, *m_avSupport);
    delete cone;
  }

  G4RotationMatrix rot;
  G4ThreeVector place;
  place.setZ(ServiceOffset * cm);
  G4Transform3D transform(rot, place);
  m_avSupport->MakeImprint(lv, transform, 0, m_overlapCheck);

  m_avBarrelCable = buildBarrelCable();
  G4ThreeVector placeBarrelCable;
  for (unsigned int i = 0; i < totStaves; ++i)
  {
    float phi = (2.0 * M_PI / totStaves) * i;
    placeBarrelCable.setX((BarrelRadius - 1) * std::cos(phi) * cm);
    placeBarrelCable.setY((BarrelRadius - 1) * std::sin(phi) * cm);
    //placeBarrelCable.setZ((-1*(BarrelLength/2 - ServiceOffset))*cm);
    G4RotationMatrix rotBarrelCable;
    rotBarrelCable.rotateZ(phi + (-90. * degToRad));
    G4Transform3D transformBarrelCable(rotBarrelCable, placeBarrelCable);
    m_avBarrelCable->MakeImprint(lv, transformBarrelCable, 0, m_overlapCheck);
  }

  m_avL0Cable = buildL0Cable();
  m_avL1Cable = buildL1Cable();
  m_avL2Cable = buildL2Cable();

  for (unsigned int iLayer = 0; iLayer < PHG4MvtxDefs::kNLayers; ++iLayer)
  {
    for (unsigned int iStave = 0; iStave < nStaves[iLayer]; ++iStave)
    {
      G4RotationMatrix rotCable;
      G4ThreeVector placeCable;
      float phi = (2.0 * M_PI / nStaves[iLayer]) * iStave;
      placeCable.setX(std::cos(phi) * cm);
      placeCable.setY(std::sin(phi) * cm);
      placeCable.setZ((ServiceOffset) *cm);
      rotCable.rotateZ(phi + ((-90. + cableRotate[iLayer]) * degToRad));
      G4Transform3D transformCable(rotCable, placeCable);
      if (iLayer == 0) m_avL0Cable->MakeImprint(lv, transformCable, 0, m_overlapCheck);
      if (iLayer == 1) m_avL1Cable->MakeImprint(lv, transformCable, 0, m_overlapCheck);
      if (iLayer == 2) m_avL2Cable->MakeImprint(lv, transformCable, 0, m_overlapCheck);
    }
  }
}
