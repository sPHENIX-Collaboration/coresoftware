#include "PHG4OuterTrackerDetector.h"

#include "PHG4OuterTrackerDefs.h"
#include "PHG4OuterTrackerDisplayAction.h"
#include "PHG4OuterTrackerSubsystem.h"

#include <outertracker/CylinderGeomOuterTracker.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Detector.h>                    // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>               // for PHG4DisplayAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>                   // for PHNodeIterator
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>            // for G4GDMLReadStructure
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4RotationMatrix.hh>               // for G4RotationMatrix
#include <Geant4/G4String.hh>                       // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>                  // for G4ThreeVector
#include <Geant4/G4Types.hh>                        // for G4double
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4VPhysicalVolume.hh>              // for G4VPhysicalVolume

#include <boost/format.hpp>

#include <cmath>
#include <cstdio>                                  // for sprintf
#include <iostream>                                 // for operator<<, basic...
#include <memory>
#include <sstream>
#include <utility>                                  // for pair, make_pair
#include <vector>                                   // for vector, vector<>:...

using namespace std;

PHG4OuterTrackerDetector::PHG4OuterTrackerDetector(PHG4OuterTrackerSubsystem* subsys, const int layer_in, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam)
  : PHG4Detector(Node, dnam)
    , m_DisplayAction(dynamic_cast<PHG4OuterTrackerDisplayAction*>(subsys->GetDisplayAction()))
{
if (Verbosity() > 0)
  cout << "PHG4OuterTrackerDetector constructor called" << endl;

if (Verbosity() > 10)
  cout << " cm " << cm << " mm " << mm << endl;

layer = layer_in;
params = _paramsContainer->GetParameters(layer);
}

//_______________________________________________________________
//_______________________________________________________________

int PHG4OuterTrackerDetector::IsInOuterTracker(G4VPhysicalVolume* volume) const
{

  if(params->get_int_param("active"))
  {
    if (activevols.find(volume) != activevols.end())
    {
      return 1;
    }
  }
  if(params->get_int_param("absorberactive"))
  {
    if (absorbervols.find(volume) != absorbervols.end())
    {
      return -1;
    }
  }
  return 0;


  if(params->get_int_param("active"))
  {
      return 1;
  }
  if(params->get_int_param("absorberactive"))
  {
      return -1;
  }
  return 0;
}

int PHG4OuterTrackerDetector::IsBlackHole(G4VPhysicalVolume* volume) const
{
  if(params->get_int_param("blackhole"))
  {
      return 1;
  }
  return 0;
}

void PHG4OuterTrackerDetector::Construct(G4LogicalVolume* logicWorld)
{
// This is called from PHG4PhenixDetector::Construct()

  if (Verbosity() > 0)
    cout << endl
         << "PHG4OuterTrackerDetector::Construct called for OuterTracker " << endl;

  ConstructOuterTracker(logicWorld);

  // This object provides the strip center locations when given the ladder segment and strip internal cordinates in the sensor
  AddGeometryNode();
  return;
}

int PHG4OuterTrackerDetector::ConstructOuterTracker(G4LogicalVolume* ot_envelope)
{
if (Verbosity() > 0)
  {
    cout << " PHG4OuterTrackerDetector::ConstructOuterTracker:" << endl;
    cout << endl;
  }

length = params->get_double_param("ot_length");
inner_radius = params->get_double_param("ot_inner_radius");
outer_radius = params->get_double_param("ot_outer_radius");
nseg_phi = params->get_int_param("ot_nseg_phi");
nseg_z = params->get_int_param("ot_nseg_z");

  // We make a simple cylinder of silicon
  
 G4VSolid *ot_si = new G4Tubs("ot_si", inner_radius * cm, outer_radius * cm, length * cm, 0.0, 2.0 * M_PI);
 
 // needs a unique name - add layer 

  G4LogicalVolume *ot_si_logic = new G4LogicalVolume(ot_si,
						     G4Material::GetMaterial("G4_Si"),
						     (boost::format("outertracker_si_volume_%d") % layer).str()  );
  
  m_DisplayAction->AddVolume(ot_si_logic, "OuterTrackerSi");
 
G4VPhysicalVolume *ot_si_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
						    ot_si_logic,
						  (boost::format("outertracker_si_%d") % layer).str(),
						    ot_envelope, false, 0, OverlapCheck() );  

 activevols.insert(ot_si_phys);

cout << ot_si_phys->GetCopyNo() << endl;
  return 0;
}

void PHG4OuterTrackerDetector::SetDisplayProperty(G4AssemblyVolume* av)
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

void PHG4OuterTrackerDetector::SetDisplayProperty(G4LogicalVolume* lv)
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

void PHG4OuterTrackerDetector::AddGeometryNode()
{
  if(params->get_int_param("active"))
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
    // here in the detector class we have internal units, convert to cm
    // before putting into the geom object
    
    CylinderGeomOuterTracker* mygeom = new CylinderGeomOuterTracker(
								      layer,
								      inner_radius,
								      outer_radius,
								      length,
								      nseg_phi,
								      nseg_z);
									
      geo->AddLayerGeom(layer, mygeom);  // only one layer per instance

    if (Verbosity())
      {
	geo->identify();
      }
  }  //is active
}  // AddGeometryNode

