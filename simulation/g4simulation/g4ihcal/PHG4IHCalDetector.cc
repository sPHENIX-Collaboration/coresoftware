#include "PHG4IHCalDetector.h"

#include "PHG4IHCalDisplayAction.h"

#include <g4detectors/PHG4HcalDefs.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>

#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VPhysicalVolume.hh>  // for G4VPhysicalVolume
#include <Geant4/G4VSolid.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#pragma GCC diagnostic pop

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>       // for unique_ptr
#include <type_traits>  // for __decay_and_strip<>::_...
#include <vector>       // for vector, vector<>::iter...

class G4Material;
class PHCompositeNode;

PHG4IHCalDetector::PHG4IHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4IHCalDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parameters)
  , m_InnerRadius(m_Params->get_double_param("inner_radius") * cm)
  , m_OuterRadius(m_Params->get_double_param("outer_radius") * cm)
  , m_SizeZ(m_Params->get_double_param("size_z") * cm)
  , m_NumScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_Active(m_Params->get_int_param("active"))
  , m_AbsorberActive(m_Params->get_int_param("absorberactive"))
  , m_GDMPath(m_Params->get_string_param("GDMPath"))
{
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
    if (m_SteelAbsorberLogVolSet.find(volume->GetLogicalVolume()) != m_SteelAbsorberLogVolSet.end())
    {
      return -1;
    }
  }
  if (m_Active)
  {
    if (m_ScintiTileLogVolSet.find(volume->GetLogicalVolume()) != m_ScintiTileLogVolSet.end())
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
  G4Material *worldmat = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid *hcal_envelope_cylinder = new G4Tubs("IHCal_envelope_solid", m_InnerRadius, m_OuterRadius, m_SizeZ / 2., 0, 2 * M_PI);
  m_VolumeEnvelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, worldmat, "Hcal_envelope", 0, 0, 0);

  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4VPhysicalVolume *mothervol = new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)), hcal_envelope_log, "IHCalEnvelope", logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyTopVolume(mothervol);
  ConstructIHCal(hcal_envelope_log);
  return;
}

int PHG4IHCalDetector::ConstructAbsorber(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope)
{
  std::vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
  for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++)
  {
    m_DisplayAction->AddSteelVolume((*it)->GetLogicalVolume());
    m_SteelAbsorberLogVolSet.insert((*it)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it));
    m_VolumeSteel += (*it)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    ++it;
  }

  return 0;
}

int PHG4IHCalDetector::ConstructScinTiles(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope)
{
  std::vector<G4VPhysicalVolume *>::iterator it = avol->GetVolumesIterator();
  for (unsigned int i = 0; i < avol->TotalImprintedVolumes(); i++)
  {
    m_DisplayAction->AddScintiVolume((*it)->GetLogicalVolume());
    m_ScintiTileLogVolSet.insert((*it)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it));
    m_ScintiTilePhysVolMap.insert(std::make_pair(*it, ExtractLayerTowerId(*it)));
    m_VolumeScintillator += (*it)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    ++it;
  }

  return 0;
}

int PHG4IHCalDetector::ConstructIHCal(G4LogicalVolume *hcalenvelope)
{
  // import the staves from the geometry file
  std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, false);

  G4AssemblyVolume *abs_asym = reader->GetAssembly("InnerSector");      //absorber
  m_ScintiMotherAssembly = reader->GetAssembly("InnerTileAssembly90");  // scintillator

  ConstructAbsorber(abs_asym, hcalenvelope);
  ConstructScinTiles(m_ScintiMotherAssembly, hcalenvelope);

  return 0;
}

int PHG4IHCalDetector::ConsistencyCheck() const
{
  // just make sure the parameters make a bit of sense
  if (m_InnerRadius >= m_OuterRadius)
  {
    std::cout << PHWHERE << ": Inner Radius " << m_InnerRadius / cm
              << " cm larger than Outer Radius " << m_OuterRadius / cm
              << " cm" << std::endl;
    gSystem->Exit(1);
  }
  return 0;
}

void PHG4IHCalDetector::Print(const std::string &what) const
{
  std::cout << "Inner Hcal Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Volume Envelope: " << m_VolumeEnvelope / cm3 << " cm^3" << std::endl;
    std::cout << "Volume Steel: " << m_VolumeSteel / cm3 << " cm^3" << std::endl;
    std::cout << "Volume Scintillator: " << m_VolumeScintillator / cm3 << " cm^3" << std::endl;
    std::cout << "Volume Air: " << (m_VolumeEnvelope - m_VolumeSteel - m_VolumeScintillator) / cm3 << " cm^3" << std::endl;
  }
  return;
}

std::pair<int, int> PHG4IHCalDetector::GetLayerTowerId(G4VPhysicalVolume *volume) const
{
  auto it = m_ScintiTilePhysVolMap.find(volume);
  if (it != m_ScintiTilePhysVolMap.end())
  {
    return it->second;
  }
  std::cout << "could not locate volume " << volume->GetName()
            << " in Inner Hcal scintillator map" << std::endl;
  gSystem->Exit(1);
  // that's dumb but code checkers do not know that gSystem->Exit()
  // terminates, so using the standard exit() makes them happy
  exit(1);
}

std::pair<int, int> PHG4IHCalDetector::ExtractLayerTowerId(G4VPhysicalVolume *volume)
{
  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char>> tok(volume->GetName(), sep);
  boost::tokenizer<boost::char_separator<char>>::const_iterator tokeniter;
  int layer_id = -1, tower_id = -1;
  for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
  {
    if (*tokeniter == "impr")
    {
      ++tokeniter;
      if (tokeniter != tok.end())
      {
        layer_id = boost::lexical_cast<int>(*tokeniter) / 2;
        layer_id--;
        if (layer_id < 0 || layer_id >= m_NumScintiPlates)
        {
          std::cout << "invalid scintillator row " << layer_id
                    << ", valid range 0 < row < " << m_NumScintiPlates << std::endl;
          gSystem->Exit(1);
        }
      }
      else
      {
        std::cout << PHWHERE << " Error parsing " << volume->GetName()
                  << " for mother volume number " << std::endl;
        gSystem->Exit(1);
      }
      break;
    }
  }
  for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
  {
    if (*tokeniter == "pv")
    {
      ++tokeniter;
      if (tokeniter != tok.end())
      {
        tower_id = boost::lexical_cast<int>(*tokeniter);
      }
    }
  }
  return std::make_pair(layer_id, tower_id);
}
