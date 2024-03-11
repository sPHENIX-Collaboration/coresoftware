#include "PHG4OHCalDetector.h"

#include "PHG4OHCalDisplayAction.h"
#include "PHG4OHCalFieldSetup.h"

#include <g4detectors/PHG4HcalDefs.h>

#include <phparameter/PHParameters.h>

#include <phfield/PHFieldConfig.h>
#include <phfield/PHFieldUtility.h>

#include <g4gdml/PHG4GDMLConfig.hh>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <calobase/RawTowerDefs.h>           // for convert_name_...
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomC...
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>

#include <g4main/PHG4Detector.h>
#include <g4main/PHG4DisplayAction.h>
#include <g4main/PHG4Subsystem.h>
#include <g4main/PHG4Utils.h>

#include <ffamodules/CDBInterface.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4IonisParamMat.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4LogicalVolumeStore.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4MaterialTable.hh>
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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <list>
#include <memory>   // for unique_ptr
#include <utility>  // for pair, make_pair
#include <vector>   // for vector, vector<>::iter...

PHG4OHCalDetector::PHG4OHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parames, const std::string &dnam)
  : PHG4Detector(subsys, Node, dnam)
  , m_DisplayAction(dynamic_cast<PHG4OHCalDisplayAction *>(subsys->GetDisplayAction()))
  , m_Params(parames)
  , m_InnerRadius(m_Params->get_double_param("inner_radius") * cm)
  , m_OuterRadius(m_Params->get_double_param("outer_radius") * cm)
  , m_SizeZ(m_Params->get_double_param("size_z") * cm)
  , m_NumScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_GDMPath(m_Params->get_string_param("GDMPath"))
{
  gdml_config = PHG4GDMLUtility::GetOrMakeConfigNode(Node);
  assert(gdml_config);
  // changes in the parameters have to be made here
  // otherwise they will not be propagated to the node tree
  if (std::filesystem::path(m_GDMPath).extension() != ".gdml")
  {
    m_GDMPath = CDBInterface::instance()->getUrl(m_GDMPath);
    m_Params->set_string_param("GDMPath", m_GDMPath);
  }
  std::string ironfieldmap = m_Params->get_string_param("IronFieldMapPath");
  if (std::filesystem::path(ironfieldmap).extension() != ".root")
  {
    ironfieldmap = CDBInterface::instance()->getUrl(ironfieldmap);
    m_Params->set_string_param("IronFieldMapPath", ironfieldmap);
  }
  PHFieldConfig *fieldconf = findNode::getClass<PHFieldConfig>(Node, PHFieldUtility::GetDSTConfigNodeName());
  if (fieldconf->get_field_config() != PHFieldConfig::kFieldUniform)
  {
    m_FieldSetup =
        new PHG4OHCalFieldSetup(
            ironfieldmap, m_Params->get_double_param("IronFieldMapScale"),
            m_InnerRadius - 10 * cm,  // subtract 10 cm to make sure fieldmap with 2x2 grid covers it
            m_OuterRadius + 10 * cm,  // add 10 cm to make sure fieldmap with 2x2 grid covers it
            m_SizeZ / 2. + 10 * cm    // div by 2 bc G4 convention
        );
  }
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
  G4LogicalVolume *hcal_envelope_log = new G4LogicalVolume(hcal_envelope_cylinder, Air, G4String("OHCal_envelope"), nullptr, nullptr, nullptr);
  G4RotationMatrix hcal_rotm;
  hcal_rotm.rotateX(m_Params->get_double_param("rot_x") * deg);
  hcal_rotm.rotateY(m_Params->get_double_param("rot_y") * deg);
  hcal_rotm.rotateZ(m_Params->get_double_param("rot_z") * deg);
  G4VPhysicalVolume *mothervol = new G4PVPlacement(G4Transform3D(hcal_rotm, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm)), hcal_envelope_log, "OHCal", logicWorld, false, false, OverlapCheck());
  m_DisplayAction->SetMyTopVolume(mothervol);
  ConstructOHCal(hcal_envelope_log);

  // allow installing new G4 subsystem installed inside the HCal envelope via macros, in particular its support rings.
  PHG4Subsystem *mysys = GetMySubsystem();
  if (mysys)
  {
    mysys->SetLogicalVolume(hcal_envelope_log);
  }
  // disable GDML export for HCal geometries for memory saving and compatibility issues
  assert(gdml_config);
  gdml_config->exclude_physical_vol(mothervol);
  gdml_config->exclude_logical_vol(hcal_envelope_log);

  const G4MaterialTable *mtable = G4Material::GetMaterialTable();
  int nMaterials = G4Material::GetNumberOfMaterials();
  for (auto i = 0; i < nMaterials; ++i)
  {
    const G4Material *mat = (*mtable)[i];
    if (mat->GetName() == "Uniplast_scintillator")
    {
      if ((mat->GetIonisation()->GetBirksConstant()) == 0)
      {
        mat->GetIonisation()->SetBirksConstant(m_Params->get_double_param("Birk_const"));
      }
    }
  }
  if (!m_Params->get_int_param("saveg4hit"))
  {
    AddGeometryNode();
  }
  return;
}

int PHG4OHCalDetector::ConstructOHCal(G4LogicalVolume *hcalenvelope)
{
  // import the staves from the gemetry file
  std::unique_ptr<G4GDMLReadStructure> reader(new G4GDMLReadStructure());
  G4GDMLParser gdmlParser(reader.get());
  gdmlParser.SetOverlapCheck(OverlapCheck());
  if (!std::filesystem::exists(m_GDMPath))
  {
    std::cout << PHWHERE << " Outer HCal gdml file " << m_GDMPath << " not found" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  gdmlParser.Read(m_GDMPath, false);

  G4AssemblyVolume *abs_asym = reader->GetAssembly("sector");         // absorber
  m_ScintiMotherAssembly = reader->GetAssembly("tileAssembly24_90");  // tiles

  // this loop is inefficient but the assignment of the scintillator id's is much simpler when having the hcal sector
  std::vector<G4VPhysicalVolume *>::iterator it1 = abs_asym->GetVolumesIterator();
  for (unsigned int isector = 0; isector < abs_asym->TotalImprintedVolumes(); isector++)
  {
    m_DisplayAction->AddSteelVolume((*it1)->GetLogicalVolume());
    m_SteelAbsorberLogVolSet.insert((*it1)->GetLogicalVolume());
    hcalenvelope->AddDaughter((*it1));
    m_VolumeSteel += (*it1)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    std::vector<G4VPhysicalVolume *>::iterator it3 = m_ScintiMotherAssembly->GetVolumesIterator();
    unsigned int ncnt = 24 * 5 * 2;
    unsigned int ioff = isector * ncnt;
    // ok we always have to skip to the scintillators we want to add for every hcal sector
    for (unsigned int j = 0; j < ioff; j++)
    {
      ++it3;
    }
    for (unsigned int j = ioff; j < ioff + ncnt; j++)
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
  // Chimney assemblies
  G4AssemblyVolume *chimAbs_asym = reader->GetAssembly("sectorChimney");         // absorber
  m_ChimScintiMotherAssembly = reader->GetAssembly("tileAssembly24chimney_90");  // chimney tiles

  std::vector<G4VPhysicalVolume *>::iterator it2 = chimAbs_asym->GetVolumesIterator();
  //	order sector 30,31,29
  std::map<unsigned int, unsigned int> sectormap;
  sectormap.insert(std::make_pair(0, 30));
  sectormap.insert(std::make_pair(1, 31));
  sectormap.insert(std::make_pair(2, 29));
  for (unsigned int isector = 0; isector < chimAbs_asym->TotalImprintedVolumes(); isector++)
  {
    m_DisplayAction->AddChimSteelVolume((*it2)->GetLogicalVolume());
    m_SteelAbsorberLogVolSet.insert((*it2)->GetLogicalVolume());

    hcalenvelope->AddDaughter((*it2));
    m_VolumeSteel += (*it2)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    std::vector<G4VPhysicalVolume *>::iterator it4 = m_ChimScintiMotherAssembly->GetVolumesIterator();
    unsigned int ncnt = 24 * 5 * 2;
    unsigned int ioff = isector * ncnt;
    // ok we always have to skip to the scintillators we want to add for every hcal sector
    for (unsigned int j = 0; j < ioff; j++)
    {
      ++it4;
    }
    for (unsigned int j = ioff; j < ioff + ncnt; j++)
    {
      m_DisplayAction->AddScintiVolume((*it4)->GetLogicalVolume());
      m_ScintiTileLogVolSet.insert((*it4)->GetLogicalVolume());
      hcalenvelope->AddDaughter((*it4));
      m_ScintiTilePhysVolMap.insert(std::make_pair(*it4, ExtractLayerTowerId(sectormap[isector], *it4)));  // chimney sectors 29-31
      m_VolumeScintillator += (*it4)->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      ++it4;
    }
    ++it2;
  }
  //Inner HCal support ring (only the part in Outer HCal volume)
  // it only exists in the new gdml file, this check keeps the old file
// without the inner hcal support readable
  std::list<std::string> volumelist = {"RingSupport_steel_1", "RingSupport_steel_2","RingSupport_alum_1","RingSupport_alum_2","HCalRing_alum_1","HCalRing_alum_2"};
  for (auto const &volumename: volumelist)
  {
    G4LogicalVolume *logvolptr = G4LogicalVolumeStore::GetInstance()->GetVolume(volumename,false);
    if (logvolptr)
    {
      m_DisplayAction->AddSupportRingVolume(logvolptr);
      m_SteelAbsorberLogVolSet.insert(logvolptr);
      hcalenvelope->AddDaughter(reader->GetWorldVolume(volumename));
    }
  }
  for (auto &logical_vol : m_SteelAbsorberLogVolSet)
  {
    if (m_FieldSetup)  // only if we have a field defined for the steel absorber
    {
      logical_vol->SetFieldManager(m_FieldSetup->get_Field_Manager_Iron(), true);

      if (m_Params->get_int_param("field_check"))
      {
        std::cout << __PRETTY_FUNCTION__ << " : setup Field_Manager_Iron for LV "
                  << logical_vol->GetName() << " w/ # of daughter " << logical_vol->GetNoDaughters() << std::endl;
      }
    }
  }

  return 0;
}

void PHG4OHCalDetector::Print(const std::string &what) const
{
  std::cout << "Outer Hcal Detector:" << std::endl;
  if (what == "ALL" || what == "VOLUME")
  {
    std::cout << "Volume Envelope: " << m_VolumeEnvelope / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Steel: " << m_VolumeSteel / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Scintillator: " << m_VolumeScintillator / cm / cm / cm << " cm^3" << std::endl;
    std::cout << "Volume Air: " << (m_VolumeEnvelope - m_VolumeSteel - m_VolumeScintillator) / cm / cm / cm << " cm^3" << std::endl;
  }
  std::cout << "******\tm_GDMPath : " << m_GDMPath << std::endl;

  return;
}

std::tuple<int, int, int> PHG4OHCalDetector::GetRowColumnId(G4VPhysicalVolume *volume) const
{
  auto it = m_ScintiTilePhysVolMap.find(volume);
  if (it != m_ScintiTilePhysVolMap.end())
  {
    return it->second;
  }
  std::cout << "could not locate volume " << volume->GetName()
            << " in Outer Hcal scintillator map" << std::endl;
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
  return std::make_tuple(isector, row, column);
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
    itwr = 12 + itmp;
  }
  else
  {
    itwr = 11 - itmp;
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
  int rowid = -1;
  if (layer_id < 225)
  {
    rowid = layer_id + 95;
  }
  else /* if (layer_id >= 225) */
  {
    rowid = layer_id - 225;
  }
  if (isector == 29)
  {
    rowid = 45 + layer_id;
  }
  else if (isector >= 30)
  {
    rowid = 75 + layer_id;
  }
  // shift the row index up by 5
  rowid += 5;
  if (rowid > 319)
  {
    rowid -= 320;
  }
  if (rowid < 0 || rowid > 319)
  {
    std::cout << "bad rowid " << rowid << " for sector " << isector << ", layer_id " << layer_id << std::endl;
    gSystem->Exit(1);
  }

  return rowid;
}

// This is dulplicated code, we can get rid of it when we have the code to make towergeom for real data reco.
void PHG4OHCalDetector::AddGeometryNode()
{
  PHNodeIterator iter(topNode());
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_SuperDetector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_SuperDetector);
    runNode->addNode(RunDetNode);
  }
  m_TowerGeomNodeName = "TOWERGEOM_" + m_SuperDetector;
  m_RawTowerGeom = findNode::getClass<RawTowerGeomContainer>(topNode(), m_TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    m_RawTowerGeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_SuperDetector));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeom, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  double innerrad = m_Params->get_double_param(PHG4HcalDefs::innerrad);
  double thickness = m_Params->get_double_param(PHG4HcalDefs::outerrad) - innerrad;
  m_RawTowerGeom->set_radius(innerrad);
  m_RawTowerGeom->set_thickness(thickness);
  m_RawTowerGeom->set_phibins(m_Params->get_int_param(PHG4HcalDefs::n_towers));
  m_RawTowerGeom->set_etabins(m_Params->get_int_param("etabins"));
  double geom_ref_radius = innerrad + thickness / 2.;
  double phistart = m_Params->get_double_param("phistart");
  if (!std::isfinite(phistart))
  {
    std::cout << PHWHERE << " phistart is not finite: " << phistart
              << ", exiting now (this will crash anyway)" << std::endl;
    gSystem->Exit(1);
  }
  for (int i = 0; i < m_Params->get_int_param(PHG4HcalDefs::n_towers); i++)
  {
    double phiend = phistart + 2. * M_PI / m_Params->get_int_param(PHG4HcalDefs::n_towers);
    std::pair<double, double> range = std::make_pair(phiend, phistart);
    phistart = phiend;
    int tempi = i + 1;
    if (tempi >= m_Params->get_int_param(PHG4HcalDefs::n_towers))
    {
      tempi -= m_Params->get_int_param(PHG4HcalDefs::n_towers);
    }
    m_RawTowerGeom->set_phibounds(tempi, range);
  }
  double etalowbound = -m_Params->get_double_param("scinti_eta_coverage_neg");
  for (int i = 0; i < m_Params->get_int_param("etabins"); i++)
  {
    // double etahibound = etalowbound + 2.2 / get_int_param("etabins");
    double etahibound = etalowbound +
                        (m_Params->get_double_param("scinti_eta_coverage_neg") + m_Params->get_double_param("scinti_eta_coverage_pos")) / m_Params->get_int_param("etabins");
    std::pair<double, double> range = std::make_pair(etalowbound, etahibound);
    m_RawTowerGeom->set_etabounds(i, range);
    etalowbound = etahibound;
  }
  for (int iphi = 0; iphi < m_RawTowerGeom->get_phibins(); iphi++)
  {
    for (int ieta = 0; ieta < m_RawTowerGeom->get_etabins(); ieta++)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(m_SuperDetector), ieta, iphi);

      const double x(geom_ref_radius * cos(m_RawTowerGeom->get_phicenter(iphi)));
      const double y(geom_ref_radius * sin(m_RawTowerGeom->get_phicenter(iphi)));
      const double z(geom_ref_radius / tan(PHG4Utils::get_theta(m_RawTowerGeom->get_etacenter(ieta))));

      RawTowerGeom *tg = m_RawTowerGeom->get_tower_geometry(key);
      if (tg)
      {
        if (Verbosity() > 0)
        {
          std::cout << "IHCalDetector::InitRun - Tower geometry " << key << " already exists" << std::endl;
        }

        if (fabs(tg->get_center_x() - x) > 1e-4)
        {
          std::cout << "IHCalDetector::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                    << std::endl;

          return;
        }
        if (fabs(tg->get_center_y() - y) > 1e-4)
        {
          std::cout << "IHCalDetector::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                    << std::endl;
          return;
        }
        if (fabs(tg->get_center_z() - z) > 1e-4)
        {
          std::cout << "IHCalDetector::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                    << std::endl;
          return;
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << "IHCalDetector::InitRun - building tower geometry " << key << "" << std::endl;
        }

        tg = new RawTowerGeomv1(key);

        tg->set_center_x(x);
        tg->set_center_y(y);
        tg->set_center_z(z);
        m_RawTowerGeom->add_tower_geometry(tg);
      }
    }
  }
  if (Verbosity() > 0)
  {
    m_RawTowerGeom->identify();
  }
}
