#include "PHG4MvtxDetector.h"

#include "PHG4MvtxDefs.h"
#include "PHG4MvtxDisplayAction.h"

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>                   // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>      // for G4Transform3D
#include <Geant4/G4Types.hh>            // for G4double
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Tubs.hh>

#include <cmath>
#include <cstdio>    // for sprintf
#include <iostream>  // for operator<<, basic...
#include <memory>
#include <sstream>
#include <utility>  // for pair, make_pair
#include <vector>   // for vector, vector<>:...

using namespace std;

namespace mvtxGeomDef
{
  double mvtx_shell_inner_radius            =  4.8  * cm;
  double skin_thickness                     =  0.01 * cm;
  double foam_core_thickness                =  0.18 * cm;
  double mvtx_shell_length                  = 50.   * cm;
  double mvtx_shell_thickness =  skin_thickness + foam_core_thickness + skin_thickness;

  double wrap_rmin = 2.1 * cm;
  double wrap_rmax = mvtx_shell_inner_radius + mvtx_shell_thickness;
  double wrap_zlen = mvtx_shell_length;
};

PHG4MvtxDetector::PHG4MvtxDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4MvtxDisplayAction*>(subsys->GetDisplayAction()))
  , m_ParamsContainer(_paramsContainer)
  , m_StaveGeometryFile(_paramsContainer->GetParameters(PHG4MvtxDefs::GLOBAL)->get_string_param("stave_geometry_file"))
  , m_EndWheelsSideS(_paramsContainer->GetParameters(PHG4MvtxDefs::GLOBAL)->get_string_param("end_wheels_sideS"))
  , m_EndWheelsSideN(_paramsContainer->GetParameters(PHG4MvtxDefs::GLOBAL)->get_string_param("end_wheels_sideN"))
{
  //  Verbosity(2);

  if (Verbosity() > 0)
    cout << "PHG4MvtxDetector constructor called" << endl;

  if (Verbosity() > 10)
    cout << " cm " << cm << " mm " << mm << endl;
  for (int ilayer = 0; ilayer < n_Layers; ++ilayer)
  {
    const PHParameters* params = m_ParamsContainer->GetParameters(ilayer);
    m_IsLayerActive[ilayer] = params->get_int_param("active");
    m_IsLayerAbsorberActive[ilayer] = params->get_int_param("absorberactive");
    m_IsBlackHole[ilayer] = params->get_int_param("blackhole");
    m_N_staves[ilayer] = params->get_int_param("N_staves");
    m_nominal_radius[ilayer] = params->get_double_param("layer_nominal_radius");
    m_nominal_phitilt[ilayer] = params->get_double_param("phitilt");
    m_nominal_phi0[ilayer] = params->get_double_param("phi0");
  }
  /*
  const PHParameters* alpide_params = m_ParamsContainer->GetParameters(PHG4MvtxDefs::ALPIDE_SEGMENTATION);
  m_PixelX = alpide_params->get_double_param("pixel_x");
  m_PixelZ = alpide_params->get_double_param("pixel_z");
  m_PixelThickness = alpide_params->get_double_param("pixel_thickness");
  */
  if (Verbosity() > 0)
  {
    cout << "PHG4MvtxDetector constructor: making Mvtx detector. " << endl;
  }
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4MvtxDetector::IsSensor(G4VPhysicalVolume* volume) const
{
  // Is this volume one of the sensors?
  // Checks if pointer matches one of our stored sensors for this layer
  if (m_SensorPV.find(volume) != m_SensorPV.end())
  {
    if (Verbosity() > 0)
    {
      cout << " -- PHG4MvtxTDetector::IsSensor --" << endl;
      cout << " volume Name : " << volume->GetName() << endl;
      cout << " -----------------------------------------" << endl;
    }
    return 1;
  }

  return 0;
}

int PHG4MvtxDetector::IsInMvtx(G4VPhysicalVolume* volume, int& layer, int& stave) const
{
  // Does this stave belong to this layer?
  // Since the Assembly volume read from GDML does not give unique pointers
  // to sensors, we need to check the stave, which is unique
  auto iter = m_StavePV.find(volume);
  if (iter != m_StavePV.end())
  {
    tie(layer, stave) = iter->second;
    if (Verbosity() > 0)
    {
      cout << " -- PHG4MvtxDetector::IsInMvtx --" << endl;
      cout << " layer: " << layer << endl;
      cout << " stave: " << stave << endl;
      cout << " volume Name : " << volume->GetName() << endl;
      cout << " stave Name  : " << iter->first->GetName() << endl;
      cout << " -----------------------------------------" << endl;
    }
    return 1;
  }

  return 0;
}

int PHG4MvtxDetector::get_layer(int index) const
{
  // Get Mvtx layer from stave index in the Mvtx
  // Mvtx stave index start from 0, i.e index = 0 for stave 0 in layer 0
  int lay = 0;
  while (!(index < m_N_staves[lay]))
  {
    index -= m_N_staves[lay];
    lay++;
  }
  return lay;
}

int PHG4MvtxDetector::get_stave(int index) const
{
  // Get stave index in the layer from stave index in the Mvtx
  int lay = 0;
  while (!(index < m_N_staves[lay]))
  {
    index -= m_N_staves[lay];
    lay++;
  }
  return index;
}

void PHG4MvtxDetector::ConstructMe(G4LogicalVolume* logicWorld)
{
  // This is called from PHG4PhenixDetector::Construct()

  if (Verbosity() > 0)
  {
    cout << endl
         << "PHG4MvtxDetector::Construct called for Mvtx " << endl;
  }

  //Create a wrapper volume
  auto tube = new G4Tubs("sol_MVTX_Wrapper", mvtxGeomDef::wrap_rmin, mvtxGeomDef::wrap_rmax,
                         mvtxGeomDef::wrap_zlen / 2.0, -M_PI, 2.0 * M_PI);
  auto world_mat = logicWorld->GetMaterial();
  auto logicMVTX = new G4LogicalVolume(tube, world_mat, "log_MVTX_Wrapper");
  new G4PVPlacement(new G4RotationMatrix(), G4ThreeVector(), logicMVTX, "MVTX_Wrapper", logicWorld, false, 0, false);

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  // this reads in the ITS stave geometry from a file and constructs the layer from it
  ConstructMvtx(logicMVTX);
  ConstructMvtxPassiveVol(logicMVTX);

  AddGeometryNode();
  return;
}

int PHG4MvtxDetector::ConstructMvtx(G4LogicalVolume* trackerenvelope)
{
  if (Verbosity() > 0)
  {
    cout << " PHG4MvtxDetector::ConstructMvtx:" << endl;
    cout << endl;
  }
  //===================================
  // Import the stave physical volume here
  //===================================

  // import the staves from the gemetry file
  std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.Read(m_StaveGeometryFile, false);

  // figure out which assembly we want
  char assemblyname[500];
  sprintf(assemblyname, "MVTXStave");

  if (Verbosity() > 0)
  {
    cout << "Geting the stave assembly named " << assemblyname << endl;
  }
  G4AssemblyVolume* av_ITSUStave = reader->GetAssembly(assemblyname);

  for (unsigned short ilayer = 0; ilayer < n_Layers; ++ilayer)
  {
    if (m_IsLayerActive[ilayer])
    {
      if (Verbosity() > 0)
      {
        cout << endl;
        cout << " Constructing Layer " << ilayer << endl;
      }
      ConstructMvtx_Layer(ilayer, av_ITSUStave, trackerenvelope);
    }
  }
  FillPVArray(av_ITSUStave);
  SetDisplayProperty(av_ITSUStave);

  return 0;
}

int PHG4MvtxDetector::ConstructMvtx_Layer(int layer, G4AssemblyVolume* av_ITSUStave, G4LogicalVolume*& trackerenvelope)
{
  //=========================================
  // Now we populate the whole layer with the staves
  //==========================================

  int N_staves = m_N_staves[layer];
  G4double layer_nominal_radius = m_nominal_radius[layer];
  G4double phitilt = m_nominal_phitilt[layer];
  G4double phi0    = m_nominal_phi0[layer]; //YCM: azimuthal offset for the first stave

  if (N_staves < 0)
  {
    // The user did not specify how many staves to use for this layer, so we calculate the best value
    // We choose a phistep that fits N_staves into this radius, but with an arclength separation of AT LEAST arcstep
    // ideally, the radius would be chosen so that numstaves = N_staves exactly, to get the closest spacing of staves possible without overlaps
    double arcstep = 12.25;
    double numstaves = 2.0 * M_PI * layer_nominal_radius / arcstep;  // this is just to print out
    N_staves = int(2.0 * M_PI * layer_nominal_radius / arcstep);     // this is the number of staves used

    // Also suggest the ideal radius for this layer
    if (Verbosity() > 0)
    {
      cout << " Calculated N_staves for layer " /*<< layer*/
           << " layer_nominal_radius " << layer_nominal_radius
           << " ITS arcstep " << arcstep
           << " circumference divided by arcstep  " << numstaves
           << " N_staves " << N_staves
           << endl;
      cout << "A radius for this layer of " << (double) N_staves * arcstep / (2.0 * M_PI) + 0.01 << " or "
           << (double) (N_staves + 1) * arcstep / (2.0 * M_PI) + 0.01 << " would produce  perfect stave spacing" << endl;
    }
  }

  G4double phistep = get_phistep(layer);  // this produces even stave spacing
  double z_location = 0.0;

  if (Verbosity() > 0)
  {
    cout << " layer " /*<< layer*/
         << " layer_nominal_radius " << layer_nominal_radius
         << " N_staves " << N_staves
         << " phistep " << phistep
         << " phitilt " << phitilt
         << " phi0    " << phi0
         << endl;
  }

  // The stave starts out at (0,0,0) and oriented so that the sensors face upward in y
  // So we need to rotate the sensor 90 degrees before placing it using phi_offset
  // note that the gdml file uses a negative phi_offset - different coord system, apparently - the following works
  double phi_offset = M_PI / 2.0;

  for (int iphi = 0; iphi < N_staves; iphi++)
  {
    // Place the ladder segment envelopes at the correct z and phi
    // This is the azimuthal angle at which we place the stave
    // phi0 is the azimuthal offset for the first stave
    G4double phi_rotation = phi0 + (double) iphi * phistep;

    G4RotationMatrix Ra;
    G4ThreeVector Ta;

    if (Verbosity() > 0)
    {
      cout << "phi_offset = " << phi_offset << " iphi " << iphi << " phi_rotation = " << phi_rotation << " phitilt " << phitilt << endl;
    }
    // It  is first rotated in phi by the azimuthal angle phi_rotation, plus the 90 degrees needed to point the face of the sensor  at the origin,  plus the tilt (if a tilt is appropriate)

    // note - if this is layer 0-2, phitilt is the additional tilt for clearance. Otherwise it is zero
    Ra.rotateZ(phi_rotation + phi_offset + phitilt);
    // Then translated as follows

    Ta.setX(layer_nominal_radius * cos(phi_rotation));
    Ta.setY(layer_nominal_radius * sin(phi_rotation));
    Ta.setZ(z_location);

    if (Verbosity() > 0)
    {
      cout << " iphi " << iphi << " phi_rotation " << phi_rotation
           << " x " << layer_nominal_radius * cos(phi_rotation)
           << " y " << layer_nominal_radius * sin(phi_rotation)
           << " z " << z_location
           << endl;
    }
    G4Transform3D Tr(Ra, Ta);

    av_ITSUStave->MakeImprint(trackerenvelope, Tr, 0, OverlapCheck());
  }

  if (Verbosity() > 0)
  {
    cout << "This layer has a total of " << N_staves << " staves" << endl;
  }
  return 0;
}

int PHG4MvtxDetector::ConstructMvtxPassiveVol(G4LogicalVolume*& lv)
{
  if (Verbosity() > 0)
  {
    cout << " PHG4MvtxDetector::ConstructMvtxServices:" << endl;
    cout << endl;
  }

  //=======================================================
  // Add an outer shell for the MVTX - moved it from INTT PHG4InttDetector.cc
  //=======================================================
  G4LogicalVolume *mvtx_shell_outer_skin_volume = GetMvtxOuterShell(lv);
  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_outer_skin_volume,
                    "mvtx_shell_outer_skin_volume", lv, false, 0, OverlapCheck());

  //===================================
  // Construct Services geometry
  //===================================
  if ( (! m_EndWheelsSideN.empty()) && (!m_EndWheelsSideS.empty()) )
  {
    // import the end_wheels from the geometry file
    std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
    G4GDMLParser gdmlParser(reader.get());

    G4RotationMatrix Ra;
    Ra.rotateY(M_PI);
    G4ThreeVector Ta;

    Ta.setZ(-30);
    G4Transform3D TrS(Ra, Ta);

    gdmlParser.Read(m_EndWheelsSideS, false);
    G4AssemblyVolume* av_EW_S = reader->GetAssembly("EndWheelsSideA");
    av_EW_S->MakeImprint(lv, TrS, 0, OverlapCheck());

    Ta.setZ(20);
    G4Transform3D TrN(Ra, Ta);

    gdmlParser.Read(m_EndWheelsSideN, false);
    G4AssemblyVolume* av_EW_N = reader->GetAssembly("EndWheelsSideC");
    av_EW_N->MakeImprint(lv, TrN, 0, OverlapCheck());
  }
  return 0;
}

G4LogicalVolume* PHG4MvtxDetector::GetMvtxOuterShell(G4LogicalVolume*& trackerenvelope)
{
  // A Rohacell foam sandwich made of 0.1 mm thick CFRP skin and 1.8 mm Rohacell 110 foam core, it has a density of 110 kg/m**3.
  //mvtx_outer_shell
  G4Tubs* mvtx_outer_shell_tube = new G4Tubs("mvtx_outer_shell",
                                              mvtxGeomDef::mvtx_shell_inner_radius,
                                              mvtxGeomDef::mvtx_shell_inner_radius + mvtxGeomDef::mvtx_shell_thickness,
                                              mvtxGeomDef::mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);

  G4LogicalVolume *mvtx_outer_shell_volume = new G4LogicalVolume(mvtx_outer_shell_tube,
                                                                 trackerenvelope->GetMaterial(),
                                                                 "mvtx_outer_shell_volume", 0, 0, 0);


  double mvtx_shell_inner_skin_inner_radius = mvtxGeomDef::mvtx_shell_inner_radius;
  double mvtx_shell_foam_core_inner_radius  = mvtx_shell_inner_skin_inner_radius + mvtxGeomDef::skin_thickness;
  double mvtx_shell_outer_skin_inner_radius = mvtx_shell_foam_core_inner_radius  + mvtxGeomDef::foam_core_thickness;

  G4Tubs *mvtx_shell_inner_skin_tube = new G4Tubs("mvtx_shell_inner_skin",
                                                  mvtx_shell_inner_skin_inner_radius,
                                                  mvtx_shell_inner_skin_inner_radius + mvtxGeomDef::skin_thickness,
                                                  mvtxGeomDef::mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);

  G4LogicalVolume *mvtx_shell_inner_skin_volume = new G4LogicalVolume(mvtx_shell_inner_skin_tube,
                                                                      G4Material::GetMaterial("CFRP_INTT"),
                                                                      "mvtx_shell_inner_skin_volume", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_inner_skin_volume,
                    "mvtx_shell_inner_skin", mvtx_outer_shell_volume, false, 0, OverlapCheck());
  //m_DisplayAction->AddVolume(mvtx_shell_inner_skin_volume, "Rail");

  G4Tubs *mvtx_shell_foam_core_tube = new G4Tubs("mvtx_shell_foam_core",
                                                 mvtx_shell_foam_core_inner_radius,
                                                 mvtx_shell_foam_core_inner_radius + mvtxGeomDef::foam_core_thickness,
                                                 mvtxGeomDef::mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);

  G4LogicalVolume *mvtx_shell_foam_core_volume = new G4LogicalVolume(mvtx_shell_foam_core_tube,
                                                                     G4Material::GetMaterial("ROHACELL_FOAM_110"),
                                                                     "mvtx_shell_foam_core_volume", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_foam_core_volume,
                    "mvtx_shell_foam_core", mvtx_outer_shell_volume, false, 0, OverlapCheck());
  //m_DisplayAction->AddVolume(mvtx_shell_foam_core_volume, "Rail");

  G4Tubs *mvtx_shell_outer_skin_tube = new G4Tubs("mvtx_shell_outer_skin",
                                                  mvtx_shell_outer_skin_inner_radius,
                                                  mvtx_shell_outer_skin_inner_radius + mvtxGeomDef::skin_thickness,
                                                  mvtxGeomDef::mvtx_shell_length / 2.0, -M_PI, 2.0 * M_PI);

  G4LogicalVolume *mvtx_shell_outer_skin_volume = new G4LogicalVolume(mvtx_shell_outer_skin_tube,
                                                                      G4Material::GetMaterial("CFRP_INTT"),
                                                                      "mvtx_shell_outer_skin_volume", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0, 0.0), mvtx_shell_outer_skin_volume,
                    "mvtx_shell_outer_skin", mvtx_outer_shell_volume, false, 0, OverlapCheck());
  //m_DisplayAction->AddVolume(mvtx_shell_outer_skin_volume, "");

  return mvtx_outer_shell_volume;
}


void PHG4MvtxDetector::SetDisplayProperty(G4AssemblyVolume* av)
{
  //  cout <<"SetDisplayProperty - G4AssemblyVolume w/ TotalImprintedVolumes "<<av->TotalImprintedVolumes()
  //   <<"/"<<av->GetImprintsCount()<<endl;

  std::vector<G4VPhysicalVolume*>::iterator it = av->GetVolumesIterator();

  int nDaughters = av->TotalImprintedVolumes();
  for (int i = 0; i < nDaughters; ++i, ++it)
  {
    //  cout <<"SetDisplayProperty - AV["<<i<<"] = "<<(*it)->GetName()<<endl;
    G4VPhysicalVolume* pv = (*it);

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
    SetDisplayProperty(worldLogical);
  }
}

void PHG4MvtxDetector::SetDisplayProperty(G4LogicalVolume* lv)
{
  string material_name(lv->GetMaterial()->GetName());

  if (Verbosity() >= 50)
  {
    cout << "SetDisplayProperty - LV " << lv->GetName() << " built with "
         << material_name << endl;
  }
  vector<string> matname = {"SI", "KAPTON", "ALUMINUM", "Carbon", "M60J3K", "WATER"};
  bool found = false;
  for (string nam : matname)
  {
    if (material_name.find(nam) != std::string::npos)
    {
      m_DisplayAction->AddVolume(lv, nam);
      if (Verbosity() >= 50)
      {
        cout << "SetDisplayProperty - LV " << lv->GetName() << " display with " << nam << endl;
      }
      found = true;
      break;
    }
  }
  if (!found)
  {
    m_DisplayAction->AddVolume(lv, "ANYTHING_ELSE");
  }
  int nDaughters = lv->GetNoDaughters();
  for (int i = 0; i < nDaughters; ++i)
  {
    G4VPhysicalVolume* pv = lv->GetDaughter(i);

    // cout <<"SetDisplayProperty - PV["<<i<<"] = "<<pv->GetName()<<endl;

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
    SetDisplayProperty(worldLogical);
  }
}

void PHG4MvtxDetector::AddGeometryNode()
{
  int active = 0;
  for (auto& isAct : m_IsLayerActive)
  {
    active |= isAct;
  }
  if (active) // At least one layer is active
  {
    ostringstream geonode;
    geonode << "CYLINDERGEOM_" << ((m_SuperDetector != "NONE") ? m_SuperDetector : m_Detector);
    PHG4CylinderGeomContainer* geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str().c_str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units(mm), convert to cm
    // before putting into the geom object
    for (unsigned short ilayer = 0; ilayer < n_Layers; ++ilayer)
    {
      CylinderGeom_Mvtx* mygeom = new CylinderGeom_Mvtx(ilayer,
                                                        m_N_staves[ilayer],
                                                        m_nominal_radius[ilayer] / cm,
                                                        get_phistep(ilayer) / rad,
                                                        m_nominal_phitilt[ilayer] / rad,
                                                        m_nominal_phi0[ilayer] / rad
                                                        );
      geo->AddLayerGeom(ilayer, mygeom);
    }  // loop per layers
    if (Verbosity())
    {
      geo->identify();
    }
  }  //is active
}  // AddGeometryNode

void PHG4MvtxDetector::FillPVArray(G4AssemblyVolume* av)
{
  if (Verbosity() > 0)
  {
    cout << "-- FillPVArray --" << endl;
  }
  std::vector<G4VPhysicalVolume*>::iterator it = av->GetVolumesIterator();

  int nDaughters = av->TotalImprintedVolumes();
  for (int i = 0; i < nDaughters; ++i, ++it)
  {
    G4VPhysicalVolume* pv = (*it);

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();
    // we only care about the staves, which contain the sensors, not the structures
    if (pv->GetName().find("MVTXHalfStave_pv") != string::npos)
    {
      int layer = get_layer(m_StavePV.size());
      int stave = get_stave(m_StavePV.size());

      m_StavePV.insert(make_pair(pv, make_tuple(layer, stave)));

      if (Verbosity() > 0)
      {
        cout << "Mvtx layer id " << layer << endl;
        cout << "Stave in layer id " << stave << endl;
        cout << "Mvtx stave count " << m_StavePV.size() << endl;
        cout << "FillPVArray - AV[" << i << "] = " << (*it)->GetName() << endl;
        cout << "              LV[" << i << "] = " << worldLogical->GetName() << endl;
      }

      FindSensor(worldLogical);
    }
    else  // in case of stave structure
    {
      if (Verbosity() > 0)
      {
        cout << "FillPVArray - AV[" << i << "] = " << (*it)->GetName() << endl;
        cout << "              LV[" << i << "] = " << worldLogical->GetName() << endl;
      }
    }
  }
}

void PHG4MvtxDetector::FindSensor(G4LogicalVolume* lv)
{
  int nDaughters = lv->GetNoDaughters();
  for (int i = 0; i < nDaughters; ++i)
  {
    G4VPhysicalVolume* pv = lv->GetDaughter(i);
    if (Verbosity() > 0)
    {
      cout << "                 PV[" << i << "]: " << pv->GetName() << endl;
    }
    if (pv->GetName().find("MVTXSensor_") != string::npos)
    {
      m_SensorPV.insert(pv);

      if (Verbosity() > 0)
      {
        cout << "                      Adding Sensor Vol <" << pv->GetName() << " (" << m_SensorPV.size() << ")>" << endl;
      }
    }

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();

    if (Verbosity() > 0)
    {
      cout << "                 LV[" << i << "]: " << worldLogical->GetName() << endl;
    }
    FindSensor(worldLogical);
  }
}
