#include "PHG4OHCalDetector.h"

#include "PHG4OHCalDisplayAction.h"
#include "PHG4OHCalFieldSetup.h"

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
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VSolid.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <Geant4/G4GDMLParser.hh>
#include <Geant4/G4GDMLReadStructure.hh>  // for G4GDMLReadStructure
#pragma GCC diagnostic pop

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

class G4Material;
class PHCompositeNode;

using namespace std;

PHG4OHCalDetector::PHG4OHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parames, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4OHCalDisplayAction *>(subsys->GetDisplayAction()))
  , m_FieldSetup(nullptr)
  , m_Params(parames)
  , m_InnerRadius(m_Params->get_double_param("inner_radius") * cm)
  , m_OuterRadius(m_Params->get_double_param("outer_radius") * cm)
  , m_SizeZ(m_Params->get_double_param("size_z") * cm)
  , m_NumScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_GDMPath(m_Params->get_string_param("GDMPath"))
{
}

PHG4OHCalDetector::~PHG4OHCalDetector()
{
  delete m_ScintiMotherAssembly;
  delete m_ChimScintiMotherAssembly;
  delete m_FieldSetup;
}

//_______________________________________________________________
//_______________________________________________________________
int PHG4OHCalDetector::IsInOHCal(G4VPhysicalVolume *volume) const
{
  if (m_AbsorberActiveFlag)
  {
    if (m_SteelAbsorberLogVolSet.find(volume->GetLogicalVolume()) != m_SteelAbsorberLogVolSet.end())
    {
      return -1;
    }
  }
  if (m_ActiveFlag)
  {
    if (m_ScintiTileLogVolSet.find(volume->GetLogicalVolume()) != m_ScintiTileLogVolSet.end())
    {
      return 1;
    }
  }
  return 0;
}

void PHG4OHCalDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  recoConsts *rc = recoConsts::instance();
  G4Material *Air = GetDetectorMaterial(rc->get_StringFlag("WorldMaterial"));
  G4VSolid *hcal_envelope_cylinder = new G4Tubs("OHCal_envelope_solid", m_InnerRadius, m_OuterRadius, m_SizeZ / 2., 0, 2 * M_PI);
  m_VolumeEnvelope = hcal_envelope_cylinder->GetCubicVolume();
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("OHCal_envelope"), 0, 0, 0);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4VPhysicalVolume *mothervol = new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)), hcal_envelope_log, "OHCal", logicWorld, 0, false, OverlapCheck());
  m_DisplayAction->SetMyTopVolume(mothervol);
  ConstructOHCal(hcal_envelope_log);
  return;
}

int PHG4OHCalDetector::ConstructOHCal(G4LogicalVolume *hcalenvelope)
{
  // import the staves from the gemetry file
  unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  gdmlParser.Read(m_GDMPath, false);


/*
  G4AssemblyVolume *abs_asym = reader->GetAssembly("sector");             //absorber
  m_ScintiMotherAssembly = reader->GetAssembly("tileAssembly24_90");             //tiles

// this loop is inefficient but the assignment of the scintillator id's is much simpler when having the hcal sector
  vector<G4VPhysicalVolume *>::iterator it1 = abs_asym->GetVolumesIterator();
  for (unsigned int isector = 0; isector < abs_asym->TotalImprintedVolumes(); isector++)
  {
    m_DisplayAction->AddSteelVolume((*it1)->GetLogicalVolume());
    m_SteelAbsorberLogVolSet.insert((*it1)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it1));
    m_VolumeSteel += (*it1)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    vector<G4VPhysicalVolume *>::iterator it3 = m_ScintiMotherAssembly->GetVolumesIterator();
    unsigned int ncnt = 24*5*2;
    unsigned int ioff = isector*ncnt;
    // ok we always have to skip to the scintillators we want to add for every hcal sector
    for (unsigned int j = 0; j < ioff; j++)
    {
      ++it3;
    }
    for (unsigned int j = ioff; j < ioff+ncnt; j++)
    {
      m_DisplayAction->AddScintiVolume((*it3)->GetLogicalVolume());
      m_ScintiTileLogVolSet.insert((*it3)->GetLogicalVolume());
      hcalenvelope->AddDaughter((*it3));
      m_ScintiTilePhysVolMap.insert(std::make_pair(*it3, ExtractLayerTowerId(isector, *it3)));
      m_VolumeScintillator += (*it3)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      ++it3;
    }

    ++it1;
  }
*/

  G4AssemblyVolume *chimAbs_asym = reader->GetAssembly("sectorChimney");  //absorber
  m_ChimScintiMotherAssembly = reader->GetAssembly("tileAssembly24chimney_90");  //chimney tiles
/*
  vector<G4VPhysicalVolume *>::iterator it2 = chimAbs_asym->GetVolumesIterator();
  for (unsigned int isector = 0; isector < chimAbs_asym->TotalImprintedVolumes(); isector++)
  {
    m_DisplayAction->AddChimSteelVolume((*it2)->GetLogicalVolume());
    m_SteelAbsorberLogVolSet.insert((*it2)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it2));
    m_VolumeSteel += (*it2)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    vector<G4VPhysicalVolume *>::iterator it4 = m_ChimScintiMotherAssembly->GetVolumesIterator();
    unsigned int ncnt = 24*5*2;
    unsigned int ioff = isector*ncnt;
    // ok we always have to skip to the scintillators we want to add for every hcal sector
    for (unsigned int j = 0; j < ioff; j++)
    {
      ++it4;
    }
    for (unsigned int j = ioff; j < ioff+ncnt; j++)
    {
      m_DisplayAction->AddScintiVolume((*it4)->GetLogicalVolume());
      m_ScintiTileLogVolSet.insert((*it4)->GetLogicalVolume());
      hcalenvelope->AddDaughter((*it4));
      m_ScintiTilePhysVolMap.insert(std::make_pair(*it4, ExtractLayerTowerId(isector+29, *it4))); // chimney sectors 29-31
//    std::pair<int, int> bla = ExtractLayerTowerId(*it4);
//    m_ScintiTilePhysVolMap.insert(std::make_pair(*it4,make_pair(255,std::get<2>(bla)) ));
      m_VolumeScintillator += (*it4)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      ++it4;
    }
    ++it2;
  }
  return 0;
*/
  vector<G4VPhysicalVolume *>::iterator it2 = chimAbs_asym->GetVolumesIterator();
  for (unsigned int i = 0; i < chimAbs_asym->TotalImprintedVolumes(); i++){
    m_DisplayAction->AddChimSteelVolume((*it2)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it2));
    ++it2;
  }

  m_ChimScintiMotherAssembly = reader->GetAssembly("tileAssembly24chimney_90");  //chimney tiles
  vector<G4VPhysicalVolume *>::iterator it4 = m_ChimScintiMotherAssembly->GetVolumesIterator();
  unsigned int ncnt = 0;
  unsigned int tilepersec = 24*5*2;
  int nsec = 29;

  for (unsigned int isector = 0; isector < m_ChimScintiMotherAssembly->TotalImprintedVolumes(); isector++)
  {
    if (ncnt >= tilepersec)
    {
      ncnt = 0;
      nsec++;
    }
//    tuple<int, int, int> bla = ExtractLayerTowerId(nsec, *it4);
    if (nsec == 30)
    {
    cout << "nsec: " << nsec << ", ncnt: " << ncnt << endl;
      m_DisplayAction->AddScintiVolume((*it4)->GetLogicalVolume());
    m_ScintiTileLogVolSet.insert((*it4)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it4));
    m_ScintiTilePhysVolMap.insert(std::make_pair(*it4, ExtractLayerTowerId(nsec, *it4)));
    }
//    std::pair<int, int> bla = ExtractLayerTowerId(*it4);
//    m_ScintiTilePhysVolMap.insert(std::make_pair(*it4,make_pair(255,std::get<2>(bla)) ));
    m_VolumeScintillator += (*it4)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    ncnt++;
    ++it4;
  }
  std::cout << "nsec: " << nsec << endl;
  std::cout << "total number of volumes: " << m_ChimScintiMotherAssembly->TotalImprintedVolumes() << endl;
  return 0;
}

void PHG4OHCalDetector::Print(const string &what) const
{
  cout << "Outer Hcal Detector:" << endl;
  if (what == "ALL" || what == "VOLUME")
  {
    cout << "Volume Envelope: " << m_VolumeEnvelope / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Steel: " << m_VolumeSteel / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Scintillator: " << m_VolumeScintillator / cm / cm / cm << " cm^3" << endl;
    cout << "Volume Air: " << (m_VolumeEnvelope - m_VolumeSteel - m_VolumeScintillator) / cm / cm / cm << " cm^3" << endl;
  }
  cout << "******\tm_GDMPath : " << m_GDMPath << endl;

  return;
}

std::tuple<int, int, int> PHG4OHCalDetector::GetRowColumnId(G4VPhysicalVolume *volume) const
{
  auto it = m_ScintiTilePhysVolMap.find(volume);
  if (it != m_ScintiTilePhysVolMap.end())
  {
    return it->second;
  }
  cout << "could not locate volume " << volume->GetName()
       << " in Outer Hcal scintillator map" << endl;
  gSystem->Exit(1);
  // that's dumb but code checkers do not know that gSystem->Exit()
  // terminates, so using the standard exit() makes them happy
  exit(1);
}

std::tuple<int, int, int> PHG4OHCalDetector::ExtractLayerTowerId(const unsigned int isector, G4VPhysicalVolume *volume)
{
  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char>> tok(volume->GetName(), sep);
  boost::tokenizer<boost::char_separator<char>>::const_iterator tokeniter;
  int layer_id = -1;
  int tower_id = -1;
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
  int column = map_towerid(tower_id);
  int row = map_layerid(isector, layer_id);
  cout << "name: " << volume->GetName() << ", sector: " << isector << ", layer id: " << layer_id << ", row: " << row << endl;
  return std::make_tuple(isector, row, column);
//  return std::make_pair(layer_id, column);
}

// map gdml tower ids to our columns
int PHG4OHCalDetector::map_towerid(const int tower_id)
{
  // odd id's go in one direction, even id's in the other one
  // this is a shortcut to derive the observed dependency
  // commented out after this code
  int itwr = -1;
  int itmp = tower_id / 2;
  if (tower_id % 2)
  {
    itwr = 11 - itmp;
  }
  else
  {
    itwr = 12 + itmp;
  }
  return itwr;
  // here is the mapping in long form
  // if (tower_id == 23)
  // {
  //   itwr = 0;
  // }
  // else if (tower_id == 21)
  // {
  //   itwr = 1;
  // }
  // else if (tower_id ==19 )
  // {
  //   itwr = 2;
  // }
  // else if (tower_id == 17)
  // {
  //   itwr = 3;
  // }
  // else if (tower_id == 15)
  // {
  //   itwr = 4;
  // }
  // else if (tower_id == 13)
  // {
  //   itwr = 5;
  // }
  // else if (tower_id == 11)
  // {
  //   itwr = 6;
  // }
  // else if (tower_id == 9)
  // {
  //   itwr = 7;
  // }
  // else if (tower_id == 7)
  // {
  //   itwr = 8;
  // }
  // else if (tower_id == 5)
  // {
  //   itwr = 9;
  // }
  // else if (tower_id == 3)
  // {
  //   itwr = 10;
  // }
  // else if (tower_id == 1)
  // {
  //   itwr = 11;
  // }
  // else if (tower_id == 0)
  // {
  //   itwr = 12;
  // }
  // else if (tower_id == 2)
  // {
  //   itwr = 13;
  // }
  // else if (tower_id == 4)
  // {
  //   itwr = 14;
  // }
  // else if (tower_id == 6)
  // {
  //   itwr = 15;
  // }
  // else if (tower_id == 8)
  // {
  //   itwr = 16;
  // }
  // else if (tower_id == 10)
  // {
  //   itwr = 17;
  // }
  // else if (tower_id == 12)
  // {
  //   itwr = 18;
  // }
  // else if (tower_id == 14)
  // {
  //   itwr = 19;
  // }
  // else if (tower_id == 16)
  // {
  //   itwr = 20;
  // }
  // else if (tower_id == 18)
  // {
  //   itwr = 21;
  // }
  // else if (tower_id == 20)
  // {
  //   itwr = 22;
  // }
  // else if (tower_id == 22)
  // {
  //   itwr = 23;
  // }
  // else
  // {
  //   std::cout << PHWHERE << " cannot map tower " << tower_id << std::endl;
  //   gSystem->Exit(1);
  //   exit(1);
  // }
}

int PHG4OHCalDetector::map_layerid(const unsigned int isector, const int layer_id)
{
  int tmp_layer = layer_id - 10*isector; // normalize to 0-9
  if (isector>=29)
  {
    tmp_layer = layer_id - 10*(isector-29);
    return layer_id;
  }
  int rowid = 4 - tmp_layer; // lower half
  if (rowid >= 0)
  {
    if (isector <= 6)
    {
      rowid += 10*(6-isector);
    }
    else if (isector <= 28)
    {
      rowid += 10*(38-isector);
    }
    else if (isector == 29)
    {
      rowid += 10*(7);
      rowid = -1;
    }
    else if (isector == 30)
    {
      rowid += 10*(8);
    }
    else
    {
      rowid = -1;
    }
  }
  else
  { // upper half
    rowid += 10;
    if (isector <= 5)
    {
      rowid += 10*(5-isector);
    }
    else if (isector <= 28)
    {
      rowid += 10*(37-isector);
    }
    else if (isector == 29)
    {
      rowid += 10*(7);
      rowid = -1;
    }
    else if (isector == 30)
    {
      rowid += 10*(8);
    }
    else
    {
      rowid = -1;
    }

  }
  // if (isector != 6)
  // {
  //   rowid = -1;
  // }
  // switch(isector)
  // {
  // case 6:
  //   break;
  // case 5:
  //   rowid += 10;
  //   break;
  // case 4:
  //   rowid += 20;
  //   break;
  // case 3:
  //   rowid += 30;
  //   break;
  // default:
  //   rowid = -1;
  //   break;
  // }
  if (rowid < 0 || rowid > 319)
  {
    cout << "bad rowid for sector " << isector << ", layer_id " << layer_id << endl;
    rowid = 255;
  }
  // if (layer_id <= 61)
  // {
  //   rowid = 61 - layer_id;
  // }
/*
  else if (layer_id > 60 && layer_id <= 191)
  {
    rowid = 191 - layer_id + 125;
  }
  else if (layer_id > 191)
  {
    rowid = 255 - layer_id + 61;
  }
  else
  {
    std::cout << PHWHERE << " cannot map layer " << layer_id << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
*/
  return rowid;
}
