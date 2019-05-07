#include "PHG4MVTXDetector.h"

#include "PHG4MVTXDefs.h"
#include "PHG4MVTXDisplayAction.h"
#include "PHG4MVTXSubsystem.h"

#include <mvtx/CylinderGeom_MVTX.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4ReflectionFactory.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>

#include <cmath>
#include <memory>
#include <sstream>

using namespace std;

PHG4MVTXDetector::PHG4MVTXDetector(PHG4MVTXSubsystem* subsys, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam)
  : PHG4Detector(Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4MVTXDisplayAction*>(subsys->GetDisplayAction()))
  , m_ParamsContainer(_paramsContainer)
  , stave_geometry_file(_paramsContainer->GetParameters(PHG4MVTXDefs::GLOBAL)->get_string_param("stave_geometry_file"))
{
  //  Verbosity(2);

  if (Verbosity() > 0)
    cout << "PHG4MVTXDetector constructor called" << endl;

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
  }
  const PHParameters* alpide_params = m_ParamsContainer->GetParameters(PHG4MVTXDefs::ALPIDE_SEGMENTATION);
  pixel_x = alpide_params->get_double_param("pixel_x");
  pixel_z = alpide_params->get_double_param("pixel_z");
  pixel_thickness = alpide_params->get_double_param("pixel_thickness");
  if (Verbosity() > 0)
    cout << "PHG4MVTXDetector constructor: making MVTX detector. " << endl;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4MVTXDetector::IsSensor(G4VPhysicalVolume* volume) const
{
  // Is this volume one of the sensors?
  // Checks if pointer matches one of our stored sensors for this layer
  if (m_SensorPV.find(volume) != m_SensorPV.end())
  {
    if (Verbosity() > 0)
    {
      cout << " -- PHG4MVTXTDetector::IsSensor --" << endl;
      cout << " volume Name : " << volume->GetName() << endl;
      cout << " -----------------------------------------" << endl;
    }
    return 1;
  }

  return 0;
}

int PHG4MVTXDetector::IsInMVTX(G4VPhysicalVolume* volume, int& layer, int& stave) const
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
      cout << " -- PHG4MVTXDetector::IsInMVTX --" << endl;
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

int PHG4MVTXDetector::get_layer(int index) const
{
  // Get MVTX layer from stave index in the MVTX
  // MVTX stave index start from 0, i.e index = 0 for stave 0 in layer 0
  int lay = 0;
  while (!(index < m_N_staves[lay]))
  {
    index -= m_N_staves[lay];
    lay++;
  }
  return lay;
}

int PHG4MVTXDetector::get_stave(int index) const
{
  // Get stave index in the layer from stave index in the MVTX
  int lay = 0;
  while (!(index < m_N_staves[lay]))
  {
    index -= m_N_staves[lay];
    lay++;
  }
  return index;
}

void PHG4MVTXDetector::Construct(G4LogicalVolume* logicWorld)
{
  // This is called from PHG4PhenixDetector::Construct()

  if (Verbosity() > 0)
    cout << endl
         << "PHG4MVTXDetector::Construct called for MVTX " << endl;

  // the tracking layers are placed directly in the world volume, since some layers are (touching) double layers
  // this reads in the ITS stave geometry from a file and constructs the layer from it
  ConstructMVTX(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4MVTXDetector::ConstructMVTX(G4LogicalVolume* trackerenvelope)
{
  if (Verbosity() > 0)
  {
    cout << " PHG4MVTXDetector::ConstructMVTX:" << endl;
    cout << endl;
  }
  //===================================
  // Import the stave physical volume here
  //===================================

  // import the staves from the gemetry file
  std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.Read(stave_geometry_file, false);

  // figure out which assembly we want
  char assemblyname[500];
  sprintf(assemblyname, "MVTXStave");

  if (Verbosity() > 0)
    cout << "Geting the stave assembly named " << assemblyname << endl;
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
      ConstructMVTX_Layer(ilayer, av_ITSUStave, trackerenvelope);
    }
  }
  FillPVArray(av_ITSUStave);
  SetDisplayProperty(av_ITSUStave);

  return 0;
}

int PHG4MVTXDetector::ConstructMVTX_Layer(int layer, G4AssemblyVolume* av_ITSUStave, G4LogicalVolume*& trackerenvelope)
{
  //=========================================
  // Now we populate the whole layer with the staves
  //==========================================

  int N_staves = m_N_staves[layer];
  G4double layer_nominal_radius = m_nominal_radius[layer];
  G4double phitilt = m_nominal_phitilt[layer];

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
         << " phtilt " << phitilt
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
    G4double phi_rotation = (double) iphi * phistep;

    G4RotationMatrix Ra;
    G4ThreeVector Ta;

    if (Verbosity() > 0)
      cout << "phi_offset = " << phi_offset << " iphi " << iphi << " phi_rotation = " << phi_rotation << " phitilt " << phitilt << endl;

    // It  is first rotated in phi by the azimuthal angle phi_rotation, plus the 90 degrees needed to point the face of the sensor  at the origin,  plus the tilt (if a tilt is appropriate)

    Ra.rotateZ(phi_rotation + phi_offset + phitilt);  // note - if this is layer 0-2, phitilt is the additional tilt for clearance. Otherwise it is zero
    // Then translated as follows

    Ta.setX(layer_nominal_radius * cos(phi_rotation));
    Ta.setY(layer_nominal_radius * sin(phi_rotation));
    Ta.setZ(z_location);

    if (Verbosity() > 0)
      cout << " iphi " << iphi << " phi_rotation " << phi_rotation
           << " x " << layer_nominal_radius * cos(phi_rotation)
           << " y " << layer_nominal_radius * sin(phi_rotation)
           << " z " << z_location
           << endl;

    G4Transform3D Tr(Ra, Ta);

    av_ITSUStave->MakeImprint(trackerenvelope, Tr, 0, OverlapCheck());
  }

  if (Verbosity() > 0)
    cout << "This layer has a total of " << N_staves << " staves" << endl;

  return 0;
}

void PHG4MVTXDetector::SetDisplayProperty(G4AssemblyVolume* av)
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

void PHG4MVTXDetector::SetDisplayProperty(G4LogicalVolume* lv)
{
  string material_name(
      lv->GetMaterial()->GetName());

  if (Verbosity() >= 5)
    cout << "SetDisplayProperty - LV " << lv->GetName() << " built with "
         << material_name << endl;
  vector<string> matname = {"SI", "KAPTON", "ALUMINUM", "Carbon", "M60J3K", "WATER"};
  bool found = false;
  for (string nam : matname)
  {
    if (material_name.find(nam) != std::string::npos)
    {
      m_DisplayAction->AddVolume(lv, nam);
      if (Verbosity() >= 5)
        cout << "SetDisplayProperty - LV " << lv->GetName() << " display with " << nam << endl;
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

void PHG4MVTXDetector::AddGeometryNode()
{
  int active = 0;
  for (auto& isAct : m_IsLayerActive)
    active |= isAct;

  if (active)
  {
    ostringstream geonode;
    geonode << "CYLINDERGEOM_" << ((superdetector != "NONE") ? superdetector : detector_type);
    PHG4CylinderGeomContainer* geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode(), geonode.str().c_str());
    if (!geo)
    {
      geo = new PHG4CylinderGeomContainer();
      PHNodeIterator iter(topNode());
      PHCompositeNode* runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject>* newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
      runNode->addNode(newNode);
    }
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    for (unsigned short ilayer = 0; ilayer < n_Layers; ++ilayer)
    {
      CylinderGeom_MVTX* mygeom = new CylinderGeom_MVTX(ilayer,
                                                        0,
                                                        m_N_staves[ilayer],
                                                        m_nominal_radius[ilayer] / cm,
                                                        get_phistep(ilayer) / rad,
                                                        m_nominal_phitilt[ilayer] / rad,
                                                        pixel_x,
                                                        pixel_z,
                                                        pixel_thickness);
      geo->AddLayerGeom(ilayer, mygeom);
    }  // loop per layers
    if (Verbosity())
    {
      geo->identify();
    }
  }  //is active
}  // AddGeometryNode

void PHG4MVTXDetector::FillPVArray(G4AssemblyVolume* av)
{
  if (Verbosity() > 0)
    cout << "-- FillPVArray --" << endl;

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
        cout << "MVTX layer id " << layer << endl;
        cout << "Stave in layer id " << stave << endl;
        cout << "MVTX stave count " << m_StavePV.size() << endl;
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

void PHG4MVTXDetector::FindSensor(G4LogicalVolume* lv)
{
  int nDaughters = lv->GetNoDaughters();
  for (int i = 0; i < nDaughters; ++i)
  {
    G4VPhysicalVolume* pv = lv->GetDaughter(i);
    if (Verbosity() > 0)
      cout << "                 PV[" << i << "]: " << pv->GetName() << endl;

    if (pv->GetName().find("MVTXSensor_") != string::npos)
    {
      m_SensorPV.insert(pv);

      if (Verbosity() > 0)
        cout << "                      Adding Sensor Vol <" << pv->GetName() << " (" << m_SensorPV.size() << ")>" << endl;
    }

    G4LogicalVolume* worldLogical = pv->GetLogicalVolume();

    if (Verbosity() > 0)
      cout << "                 LV[" << i << "]: " << worldLogical->GetName() << endl;

    FindSensor(worldLogical);
  }
}
