#include "PHG4MvtxDetector.h"

#include "PHG4MvtxDefs.h"
#include "PHG4MvtxDisplayAction.h"
#include "PHG4MvtxSupport.h"

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>      // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Transform3D.hh>  // for G4Transform3D
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4Types.hh>            // for G4double
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume

// xerces has some shadowed variables
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#pragma GCC diagnostic pop

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
  const double wrap_rmax = ( 107.7 + 0.3 ) * mm; // 300 um Marging from SB Flange
  const double wrap_rmin = 22 * mm;
  const double wrap_smallCylR = 55.00 * mm; // CYSS Ext Nose + 905 um
  const double wrap_CYSSFlgN_Z = ( 177.5 + 12.5 ) * mm;
  const double wrap_CYSSNose_Z = -245 * mm;
  const double wrap_CYSSHead_Z = - 315 * mm;
  const double wrap_SBCyl_Z = -1800 * mm; // SB Cyl (1650 mm + 15 cm Margin)
}  // namespace mvtxGeomDef

PHG4MvtxDetector::PHG4MvtxDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4MvtxDisplayAction*>(subsys->GetDisplayAction()))
  , m_ParamsContainer(_paramsContainer)
  , m_StaveGeometryFile(_paramsContainer->GetParameters(PHG4MvtxDefs::GLOBAL)->get_string_param("stave_geometry_file"))
{
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

  if (Verbosity() > 0)
  {
    cout << "PHG4MvtxDetector constructor: making Mvtx detector. " << endl;
  }
}


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

  const G4int numZPlanes = 4;
  const G4double zPlane[numZPlanes] = { mvtxGeomDef::wrap_SBCyl_Z,
                                        mvtxGeomDef::wrap_CYSSHead_Z,
                                        mvtxGeomDef::wrap_CYSSNose_Z,
                                        mvtxGeomDef::wrap_CYSSFlgN_Z};

  const G4double rInner[numZPlanes] = { mvtxGeomDef::wrap_rmin, mvtxGeomDef::wrap_rmin,
                                        mvtxGeomDef::wrap_rmin, mvtxGeomDef::wrap_rmin };

  const G4double rOuter[numZPlanes] = { mvtxGeomDef::wrap_rmax, mvtxGeomDef::wrap_rmax,
                                        mvtxGeomDef::wrap_smallCylR,
                                        mvtxGeomDef::wrap_smallCylR };

  auto mvtxWrapSol = new G4Polycone( "sol_MVTX_Wrapper", 0, 2.0 * M_PI,
                                     numZPlanes, zPlane, rInner, rOuter );

  auto world_mat = logicWorld->GetMaterial();

  auto logicMVTX = new G4LogicalVolume( mvtxWrapSol, world_mat, "log_MVTX_Wrapper" );

  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr( Ra, Ta );
  new G4PVPlacement( Tr, logicMVTX, "MVTX_Wrapper",
                     logicWorld, false, 0, false );

  // the tracking layers are placed directly in the world volume,
  // since some layers are (touching) double layers
  // this reads in the ITS stave geometry from a file and constructs the layer from it
  ConstructMvtx( logicMVTX );
  ConstructMvtxPassiveVol( logicMVTX );

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
  G4double phi0 = m_nominal_phi0[layer];  //YCM: azimuthal offset for the first stave

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

  //Now construct EWs, service barrel, CYSS, cones and cables
  PHG4MvtxSupport* mvtxSupportSystem = new PHG4MvtxSupport(m_DisplayAction, OverlapCheck());
  mvtxSupportSystem->ConstructMvtxSupport(lv);

  delete mvtxSupportSystem;

  return 0;
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
  for (const string& nam : matname)
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
  if (active)  // At least one layer is active
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
                                                        m_nominal_phi0[ilayer] / rad);
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
