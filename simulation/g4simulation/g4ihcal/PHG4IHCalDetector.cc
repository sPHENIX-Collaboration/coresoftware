#include "PHG4IHCalDetector.h"

#include "g4detectors/PHG4HcalDefs.h"
#include "PHG4IHCalDisplayAction.h"
#include "PHG4IHCalSubsystem.h"


#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VSolid.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Circular_kernel_intersections.h>
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Object.h>
#include <CGAL/point_generators_2.h>
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#include <Geant4/G4GDMLReadParamvol.hh>
#pragma GCC diagnostic pop

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>

class PHCompositeNode;

typedef CGAL::Circle_2<PHG4IHCalDetector::Circular_k> Circle_2;
typedef CGAL::Circular_arc_point_2<PHG4IHCalDetector::Circular_k> Circular_arc_point_2;
typedef CGAL::Line_2<PHG4IHCalDetector::Circular_k> Line_2;
typedef CGAL::Segment_2<PHG4IHCalDetector::Circular_k> Segment_2;

using namespace std;

// there is still a minute problem for very low tilt angles where the scintillator
// face touches the boundary instead of the corner, subtracting 1 permille from the total
// scintilator length takes care of this
//-m/s-static double subtract_from_scinti_x = 0.1 * mm;

PHG4IHCalDetector::PHG4IHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4IHCalDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_ScintiMotherAssembly(nullptr)
  , m_InnerRadius(m_Params->get_double_param("inner_radius") * cm)
  , m_OuterRadius(m_Params->get_double_param("outer_radius") * cm)
  , m_SizeZ(m_Params->get_double_param("size_z") * cm)
  , m_ScintiTileX(NAN)
  , m_ScintiTileXLower(NAN)
  , m_ScintiTileXUpper(NAN)
  , m_ScintiTileZ(m_SizeZ)
  , m_ScintiTileThickness(m_Params->get_double_param("scinti_tile_thickness") * cm)
  , m_ScintiInnerGap(m_Params->get_double_param("scinti_inner_gap") * cm)
  , m_ScintiOuterGap(m_Params->get_double_param("scinti_outer_gap") * cm)
  , m_ScintiOuterRadius(m_Params->get_double_param("scinti_outer_radius") * cm)
  , m_TiltAngle(m_Params->get_double_param("tilt_angle") * deg)
  , m_EnvelopeInnerRadius(m_InnerRadius)
  , m_EnvelopeOuterRadius(m_OuterRadius)
  , m_EnvelopeZ(m_SizeZ)
  , m_VolumeEnvelope(NAN)
  , m_VolumeSteel(NAN)
  , m_VolumeScintillator(NAN)
  , m_NumScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_NumScintiTilesPos(m_Params->get_int_param("n_scinti_tiles_pos"))
  , m_NumScintiTilesNeg(m_Params->get_int_param("n_scinti_tiles_neg"))
  , m_Active(m_Params->get_int_param("active"))
  , m_AbsorberActive(m_Params->get_int_param("absorberactive"))
  , m_Layer(0)
  , m_ScintiLogicNamePrefix("HcalInnerScinti")
  , m_GDMPath(m_Params->get_string_param("GDMPath"))
{
  // n_scinti_tiles takes precedence

  int nTiles = 2 * m_Params->get_int_param(PHG4HcalDefs::n_scinti_tiles);
  if (nTiles <= 0)
  {
    nTiles = m_NumScintiTilesPos + m_NumScintiTilesNeg;
  }
  else
  {
    m_NumScintiTilesPos = nTiles / 2;
    m_Params->set_int_param("n_scinti_tiles_pos", nTiles / 2);
    m_NumScintiTilesNeg = nTiles / 2;
    m_Params->set_int_param("n_scinti_tiles_neg", nTiles / 2);
  }

  // allocate memory for scintillator plates
  m_ScintiTilesVec.assign(nTiles, static_cast<G4VSolid *>(nullptr));
}

PHG4IHCalDetector::~PHG4IHCalDetector()
{
  delete m_ScintiMotherAssembly;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4IHCalDetector::IsInIHCal(G4VPhysicalVolume *volume) const
{
  if (m_AbsorberActive)
  {
    if (m_SteelAbsorberPhysVolSet.find(volume) != m_SteelAbsorberPhysVolSet.end())
    {
      return -1;
    }
  }
  if (m_Active)
  {
    if (m_ScintiTilePhysVolMap.find(volume) != m_ScintiTilePhysVolMap.end())
    {
      return 1;
    }
  }
  return 0;
}

// Construct the envelope and the call the
// actual inner hcal construction
void PHG4IHCalDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  recoConsts *rc = recoConsts::instance();
  G4Material *Air = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid *hcal_envelope_cylinder = new G4Tubs("IHCal_envelope_solid", m_EnvelopeInnerRadius, m_EnvelopeOuterRadius, m_EnvelopeZ / 2., 0, 2 * M_PI);
  m_VolumeEnvelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("Hcal_envelope"), 0, 0, 0);

  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4VPhysicalVolume *mothervol = new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)), hcal_envelope_log, "IHCalEnvelope", logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyTopVolume(mothervol);
  ConstructIHCal(hcal_envelope_log);

  vector<G4VPhysicalVolume *>::iterator it = m_ScintiMotherAssembly->GetVolumesIterator();
  for (unsigned int i = 0; i < m_ScintiMotherAssembly->TotalImprintedVolumes(); i++){

    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char>> tok((*it)->GetName(), sep);
    boost::tokenizer<boost::char_separator<char>>::const_iterator tokeniter;
    int layer_id = -1, tower_id = -1;
    for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter){
      if (*tokeniter == "impr"){
        ++tokeniter;
        if (tokeniter != tok.end()){
          layer_id = boost::lexical_cast<int>(*tokeniter)/2;
	  layer_id--;
          if (layer_id < 0 || layer_id >= m_NumScintiPlates){
            cout << "invalid scintillator row " << layer_id
                 << ", valid range 0 < row < " << m_NumScintiPlates << endl;
            gSystem->Exit(1);
          }
        } else {
          cout << PHWHERE << " Error parsing " << (*it)->GetName()
               << " for mother volume number " << endl;
          gSystem->Exit(1);
        }
        break;
      }
    }
    for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter){
      if (*tokeniter == "pv"){
        ++tokeniter;
        if (tokeniter != tok.end()){
          tower_id = boost::lexical_cast<int>(*tokeniter);
	}
      }
    }
    pair<int, int> layer_twr = make_pair(layer_id, tower_id);
    m_ScintiTilePhysVolMap.insert(pair<G4VPhysicalVolume *, pair<int, int>>(*it, layer_twr));

    ++it;
  }

  return;
}

int PHG4IHCalDetector::ConstructAbsorber(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope)
{

  vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
  for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++){
    m_DisplayAction->AddSteelVolume((*it)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it));
    ++it;
  }

  return 0;
}

int PHG4IHCalDetector::ConstructScinTiles(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope)
{

  vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
  for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++){
    m_DisplayAction->AddScintiVolume((*it)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it));
    ++it;
  }

  return 0;
}

int PHG4IHCalDetector::ConstructIHCal(G4LogicalVolume *hcalenvelope)
{

  // import the staves from the gemetry file
  unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, false);

  G4AssemblyVolume* abs_asym = reader->GetAssembly("InnerSector"); //absorber
  m_ScintiMotherAssembly = reader->GetAssembly("InnerTileAssembly90"); // scintillator

  ConstructAbsorber(abs_asym, hcalenvelope);
  ConstructScinTiles(m_ScintiMotherAssembly, hcalenvelope);

  return 0;
}

int PHG4IHCalDetector::ConsistencyCheck() const
{
  // just make sure the parameters make a bit of sense
  if (m_InnerRadius >= m_OuterRadius)
  {
    cout << PHWHERE << ": Inner Radius " << m_InnerRadius / cm
         << " cm larger than Outer Radius " << m_OuterRadius / cm
         << " cm" << endl;
    gSystem->Exit(1);
  }
  return 0;
}

void PHG4IHCalDetector::Print(const string &what) const
{
  cout << "Inner Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Envelope: " << m_VolumeEnvelope / cm3 << " cm^3" << endl;
    cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << endl;
    cout << "Volume Air: " << (m_VolumeEnvelope - m_VolumeSteel - m_VolumeScintillator) / cm3 << " cm^3" << endl;
  }
  cout<<"******\tm_GDMPath : "<<m_GDMPath<<endl;

  return;
}

std::pair<int, int> PHG4IHCalDetector::GetLayerTowerId(G4VPhysicalVolume *volume) const
{
  auto it = m_ScintiTilePhysVolMap.find(volume);
  if (it != m_ScintiTilePhysVolMap.end())
  {
    return it->second;
  }
  cout << "could not locate volume " << volume->GetName()
       << " in Inner Hcal scintillator map" << endl;
  gSystem->Exit(1);
  // that's dumb but code checkers do not know that gSystem->Exit()
  // terminates, so using the standard exit() makes them happy
  exit(1);
}
