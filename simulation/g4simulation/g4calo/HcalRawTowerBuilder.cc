#include "HcalRawTowerBuilder.h"

#include <calobase/RawTower.h>  // for RawTower
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>           // for convert_name_...
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomC...
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <g4detectors/PHG4HcalDefs.h>

#include <phparameter/PHParameterInterface.h>  // for PHParameterIn...
#include <phparameter/PHParameters.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cmath>    // for fabs, NAN, cos
#include <cstdlib>  // for exit
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>  // for allocator_tra...
#include <stdexcept>
#include <utility>  // for make_pair, pair

HcalRawTowerBuilder::HcalRawTowerBuilder(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

int HcalRawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  if (m_InputDetector.empty() || m_OutputDetector.empty())
  {
    std::cout << PHWHERE
              << " Detector name not set, use HcalRawTowerBuilder::Detector(string) to set, exiting"
              << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, exiting" << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  std::string paramnodename = "TOWERPARAM_" + m_OutputDetector;
  try
    {
      CreateNodes(topNode);
    }

  catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  // order first default,
  // then parameter from g4detector on node tree
  ReadParamsFromNodeTree(topNode);
  // then macro setting
  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_OutputDetector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_OutputDetector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  std::string geonodename = "TOWERGEO_" + m_OutputDetector;

  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", m_OutputDetector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(m_OutputDetector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);
  m_TowerEnergySrc = get_int_param("tower_energy_source");
  m_Emin = get_double_param("emin");
  m_NcellToTower = get_int_param("n_scinti_plates_per_tower");
  if (!m_TowerDecalFactors.empty())
  {
    SetTowerDecalFactors();
  }
  if (Verbosity() >= 1)
  {
    std::cout << "HcalRawTowerBuilder::InitRun :";
    if (m_TowerEnergySrc == kEnergyDeposition)
    {
      std::cout << "save Geant4 energy deposition in towers" << std::endl;
    }
    else if (m_TowerEnergySrc == kLightYield)
    {
      std::cout << "save light yield in towers" << std::endl;
    }
    else if (m_TowerEnergySrc == kIonizationEnergy)
    {
      std::cout << "save ionization energy in towers" << std::endl;
    }
    else if (m_TowerEnergySrc == kRawLightYield)
    {
      std::cout << "save raw (pre-Mephi map) light yield in towers" << std::endl;
    }
    else
    {
      std::cout << "unknown energy source" << std::endl;
    }
  }
  m_TowerGeomNodeName = "TOWERGEOM_" + m_OutputDetector;
  m_RawTowerGeom = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!m_RawTowerGeom)
  {
    m_RawTowerGeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(m_InputDetector));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeom, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  double innerrad = get_double_param(PHG4HcalDefs::innerrad);
  double thickness = get_double_param(PHG4HcalDefs::outerrad) - innerrad;
  m_RawTowerGeom->set_radius(innerrad);
  m_RawTowerGeom->set_thickness(thickness);
  m_RawTowerGeom->set_phibins(get_int_param(PHG4HcalDefs::n_towers));
  m_RawTowerGeom->set_etabins(get_int_param("etabins"));
  double geom_ref_radius = innerrad + thickness / 2.;
  double phistart = get_double_param("phistart");
  if (!std::isfinite(phistart))
  {
    std::cout << PHWHERE << " phistart is not finite: " << phistart
              << ", exiting now (this will crash anyway)" << std::endl;
    gSystem->Exit(1);
  }
  for (int i = 0; i < get_int_param(PHG4HcalDefs::n_towers); i++)
  {
    double phiend = phistart + 2. * M_PI / get_int_param(PHG4HcalDefs::n_towers);
    std::pair<double, double> range = std::make_pair(phiend, phistart);
    phistart = phiend;
    m_RawTowerGeom->set_phibounds(i, range);
  }
  //double etalowbound = -1.1;
  double etalowbound = -get_double_param("scinti_eta_coverage_neg");
  for (int i = 0; i < get_int_param("etabins"); i++)
  {
    //double etahibound = etalowbound + 2.2 / get_int_param("etabins");
    double etahibound = etalowbound +
                        (get_double_param("scinti_eta_coverage_neg") + get_double_param("scinti_eta_coverage_pos")) / get_int_param("etabins");
    std::pair<double, double> range = std::make_pair(etalowbound, etahibound);
    m_RawTowerGeom->set_etabounds(i, range);
    etalowbound = etahibound;
  }
  for (int iphi = 0; iphi < m_RawTowerGeom->get_phibins(); iphi++)
  {
    for (int ieta = 0; ieta < m_RawTowerGeom->get_etabins(); ieta++)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(m_InputDetector), ieta, iphi);

      const double x(geom_ref_radius * cos(m_RawTowerGeom->get_phicenter(iphi)));
      const double y(geom_ref_radius * sin(m_RawTowerGeom->get_phicenter(iphi)));
      const double z(geom_ref_radius / tan(PHG4Utils::get_theta(m_RawTowerGeom->get_etacenter(ieta))));

      RawTowerGeom *tg = m_RawTowerGeom->get_tower_geometry(key);
      if (tg)
      {
        if (Verbosity() > 0)
        {
          std::cout << "HcalRawTowerBuilder::InitRun - Tower geometry " << key << " already exists" << std::endl;
        }

        if (fabs(tg->get_center_x() - x) > 1e-4)
        {
          std::cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                    << std::endl;

          return Fun4AllReturnCodes::ABORTRUN;
        }
        if (fabs(tg->get_center_y() - y) > 1e-4)
        {
          std::cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                    << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
        if (fabs(tg->get_center_z() - z) > 1e-4)
        {
          std::cout << "HcalRawTowerBuilder::InitRun - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                    << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << "HcalRawTowerBuilder::InitRun - building tower geometry " << key << "" << std::endl;
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

  int m = m_DecalArray[0].size();
  int n = m_DecalArray.size();

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      m_DecalArray[i][j] = 1.;
    }
  }

  if (!m_DeCalibrationFileName.empty())
  {
    if (std::filesystem::exists(m_DeCalibrationFileName))
    {
      std::ifstream decalibrate_tower;
      decalibrate_tower.open(m_DeCalibrationFileName, std::ifstream::in);
      if (decalibrate_tower.is_open())
      {
        while (!decalibrate_tower.eof())
        {
          int etabin = -1;
          int phibin = -1;
          for (int i = 0; i < n; i++)
          {
            for (int j = 0; j < m; j++)
            {
              decalibrate_tower >> etabin >> phibin >> m_DecalArray[i][j];
              if (!std::isfinite(m_DecalArray[i][j]))
              {
                std::cout << "Calibration constant at etabin " << etabin
                          << ", phibin " << phibin << " in " << m_DeCalibrationFileName
                          << " is not finite: " << m_DecalArray[i][j] << std::endl;
                gSystem->Exit(1);
                exit(1);
              }
            }
          }
        }
        decalibrate_tower.close();
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int HcalRawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  /* decalibration occurs if user supplies a non empty decalMap.txt
    file, otherwise code will proceed with no de-calibration (as is)
*/

  double cell_weight = 0.0;
  if (Verbosity() > 3)
  {
    std::cout << PHWHERE << "Process event entered" << std::endl;
  }

  //load get TowerInfoContainer node from node tree:
  TowerInfoContainer*  m_TowerInfoContainer = findNode::getClass<TowerInfoContainer>(topNode,m_TowerInfoNodeName);
  if (!m_TowerInfoContainer && m_UseTowerInfo > 0)
    {
      std::cout << PHWHERE << "TowerInfoContainer Node missing, doing nothing." << std::endl;
      exit(1);
    }

  // get cells
  std::string cellnodename = "G4CELL_" + m_InputDetector;
  PHG4CellContainer *slats = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!slats)
  {
    std::cout << PHWHERE << " Node " << cellnodename
              << " missing, quitting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  // loop over all slats in an event
  PHG4CellContainer::ConstIterator cell_iter;
  PHG4CellContainer::ConstRange cell_range = slats->getCells();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
       ++cell_iter)
  {
    PHG4Cell *cell = cell_iter->second;

    short twrrow = get_tower_row(PHG4CellDefs::ScintillatorSlatBinning::get_row(cell->get_cellid()));

    if (m_TowerEnergySrc == kEnergyDeposition)
    {
      cell_weight = cell->get_edep();
    }
    else if (m_TowerEnergySrc == kLightYield)
    {
      cell_weight = cell->get_light_yield();
    }
    else if (m_TowerEnergySrc == kIonizationEnergy)
    {
      cell_weight = cell->get_eion();
    }
    else if (m_TowerEnergySrc == kRawLightYield)
    {
      cell_weight = cell->get_raw_light_yield();
    }
    else
    {
      std::cout << Name() << ": unknown tower energy source "
                << m_TowerEnergySrc << std::endl;
      gSystem->Exit(1);
      exit(1);
    }

    cell_weight *= m_DecalArray.at(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid())).at(PHG4CellDefs::ScintillatorSlatBinning::get_row(cell->get_cellid()));





    if (m_UseTowerInfo  != 1)
      {

	RawTower *tower = m_Towers->getTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow);
	if (!tower)
	  {
	    tower = new RawTowerv1();
	    tower->set_energy(0.0);
	    m_Towers->AddTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow, tower);
	  }
	
	
	tower->add_ecell(cell->get_cellid(), cell_weight);
	
	PHG4Cell::ShowerEdepConstRange range = cell->get_g4showers();
	for (PHG4Cell::ShowerEdepConstIterator shower_iter = range.first;
	     shower_iter != range.second;
	     ++shower_iter)
	  {
	    tower->add_eshower(shower_iter->first, shower_iter->second);
	  }
	tower->set_energy(tower->get_energy() + cell_weight);
    }
    
    if (m_UseTowerInfo > 0)
      {
	TowerInfo *towerinfo;
	unsigned int etabin = PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid());
	unsigned int phibin = twrrow;
	
	unsigned int towerkey = (etabin << 16U) + phibin;
	unsigned int towerindex = m_TowerInfoContainer->decode_key(towerkey);

        towerinfo = m_TowerInfoContainer->at(towerindex);
        if (!towerinfo)
        {
          std::cout << __PRETTY_FUNCTION__ << ": missing towerkey = " << towerkey << " in m_TowerInfoContainer!";
          exit(1);
        }
        else
        {
          towerinfo->set_energy(towerinfo->get_energy() + cell_weight);
        }
      }
  }

  double towerE = 0;
  if (m_ChkEnergyConservationFlag)
  {
    double cellE = slats->getTotalEdep();
    towerE = m_Towers->getTotalEdep();
    if (fabs(cellE - towerE) / cellE > 1e-5)
    {
      std::cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
                << cellE - towerE << std::endl;
    }
  }

  if (Verbosity())
  {
    towerE = m_Towers->getTotalEdep();
  }

  m_Towers->compress(m_Emin);
  if (Verbosity())
  {
    std::cout << "Energy lost by dropping towers with less than " << m_Emin
              << " energy, lost energy: " << towerE - m_Towers->getTotalEdep()
              << std::endl;
    m_Towers->identify();
    RawTowerContainer::ConstRange begin_end = m_Towers->getTowers();
    RawTowerContainer::ConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter)
    {
      iter->second->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void HcalRawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, exiting." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_OutputDetector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_OutputDetector);
    dstNode->addNode(DetNode);
  }





  if (m_UseTowerInfo != 1)
    {
      // Create the tower nodes on the tree
      if (m_SimTowerNodePrefix.empty())
	{
	  // no prefix, consistent with older convension
	  m_TowerNodeName = "TOWER_" + m_OutputDetector;
	}
      else
	{
	  m_TowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_OutputDetector;
	}
      m_Towers = findNode::getClass<RawTowerContainer>(DetNode, m_TowerNodeName);
      if (!m_Towers)
	{
	  m_Towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(m_InputDetector));
	  
	  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_Towers, m_TowerNodeName, "PHObject");
	  DetNode->addNode(towerNode);
	}
    }


  if (m_UseTowerInfo > 0 )
    {


     if (m_SimTowerNodePrefix.empty())
	{
	  // no prefix, consistent with older convension
	  m_TowerInfoNodeName = "TOWERINFO_" + m_OutputDetector;
	}
      else
	{
	  m_TowerInfoNodeName = "TOWERINFO_" + m_SimTowerNodePrefix + "_" + m_OutputDetector;
	}


     TowerInfoContainer* m_TowerInfoContainer = findNode::getClass<TowerInfoContainer>(DetNode,m_TowerInfoNodeName);

     RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(m_InputDetector);
     if (m_TowerInfoContainer == nullptr)
       {
	 TowerInfoContainerv1::DETECTOR detec;
	 if (caloid == RawTowerDefs::CalorimeterId::CEMC)
	   {
	     detec = TowerInfoContainer::DETECTOR::EMCAL;
	   }
	 else if (caloid == RawTowerDefs::CalorimeterId::HCALIN || caloid == RawTowerDefs::CalorimeterId::HCALOUT)
	   {
	     detec = TowerInfoContainer::DETECTOR::HCAL;
	   }
	 else
	   {
	     std::cout << PHWHERE << "Detector not implemented into the TowerInfoContainer object, defaulting to HCal implementation." << std::endl;
	     detec = TowerInfoContainer::DETECTOR::HCAL;
	   }
	 m_TowerInfoContainer = new TowerInfoContainerv1(detec);
	 PHIODataNode<PHObject> *TowerInfoNode = new PHIODataNode<PHObject>(m_TowerInfoContainer, m_TowerInfoNodeName, "PHObject");
	 DetNode->addNode(TowerInfoNode);
       }






    }


  return;
}

short HcalRawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow / m_NcellToTower;
  return twrrow;
}

void HcalRawTowerBuilder::SetDefaultParameters()
{
  set_default_int_param(PHG4HcalDefs::scipertwr, 5);
  set_default_int_param("tower_energy_source", kLightYield);
  set_default_int_param(PHG4HcalDefs::n_towers, 64);
  set_default_int_param("etabins", 24);

  set_default_double_param("emin", 1.e-6);
  set_default_double_param(PHG4HcalDefs::outerrad, NAN);
  set_default_double_param(PHG4HcalDefs::innerrad, NAN);

  set_default_double_param("scinti_eta_coverage_neg", 1.1);
  set_default_double_param("scinti_eta_coverage_pos", 1.1);
  set_default_double_param("phistart", NAN);
}

void HcalRawTowerBuilder::ReadParamsFromNodeTree(PHCompositeNode *topNode)
{
  PHParameters *pars = new PHParameters("temp");
  // we need the number of scintillator plates per tower
  std::string geonodename = "G4GEOPARAM_" + m_InputDetector;
  PdbParameterMapContainer *saveparams = findNode::getClass<PdbParameterMapContainer>(topNode, geonodename);
  if (!saveparams)
  {
    std::cout << "could not find " << geonodename << std::endl;
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Print("NODETREE");
    return;
  }
  pars->FillFrom(saveparams, 0);
  set_int_param(PHG4HcalDefs::scipertwr, pars->get_int_param(PHG4HcalDefs::scipertwr));
  set_int_param(PHG4HcalDefs::n_towers, pars->get_int_param(PHG4HcalDefs::n_towers));
  set_double_param(PHG4HcalDefs::innerrad, pars->get_double_param(PHG4HcalDefs::innerrad));
  set_double_param(PHG4HcalDefs::outerrad, pars->get_double_param(PHG4HcalDefs::outerrad));

  int nTiles = 2 * pars->get_int_param(PHG4HcalDefs::n_scinti_tiles);
  int nPhislices = pars->get_int_param(PHG4HcalDefs::scipertwr) * pars->get_int_param(PHG4HcalDefs::n_towers);
  if (nTiles <= 0)
  {
    nTiles = pars->get_int_param(PHG4HcalDefs::n_scinti_tiles_pos) + pars->get_int_param(PHG4HcalDefs::n_scinti_tiles_neg);
    set_double_param("scinti_eta_coverage_neg", pars->get_double_param("scinti_eta_coverage_neg"));
    set_double_param("scinti_eta_coverage_pos", pars->get_double_param("scinti_eta_coverage_pos"));
  }
  set_int_param("etabins", nTiles);
  m_DecalArray.resize(nTiles, std::vector<double>(nPhislices));

  delete pars;
  return;
}

void HcalRawTowerBuilder::set_cell_decal_factor(const int etabin, const int phibin, const double d)
{
  m_DecalArray.at(etabin).at(phibin) = d;
}

void HcalRawTowerBuilder::SetTowerDecalFactors()
{
  for (auto iter = m_TowerDecalFactors.begin(); iter != m_TowerDecalFactors.end(); ++iter)
  {
    set_tower_decal_factor_real(iter->first.first, iter->first.second, iter->second);
  }
}

void HcalRawTowerBuilder::set_tower_decal_factor(const int etabin, const int phibin, const double d)
{
  // since we do not have the number of scintillators per tower at this point
  // the decal values are cached in m_TowerDecalFactors to be set during the InitRun
  std::pair<int, int> etaphi = std::make_pair(etabin, phibin);
  m_TowerDecalFactors[etaphi] = d;
}

void HcalRawTowerBuilder::set_tower_decal_factor_real(const int etabin, const int phibin, const double d)
{
  for (int i = 0; i < m_NcellToTower; i++)
  {
    int istart = phibin * m_NcellToTower + i;
    m_DecalArray.at(etabin).at(istart) = d;
  }
}

void HcalRawTowerBuilder::Print(const std::string & /*what*/) const
{
  std::cout << Name() << std::endl;
  PHParameterInterface::Print();
}
