#include "RawTowerBuilder.h"

#include <calobase/RawTower.h>  // for RawTower
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>           // for encode_towerid
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomC...
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>





#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <boost/io/ios_state.hpp>

#include <cmath>      // for fabs, tan, atan2
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>  // for pair, make_pair

RawTowerBuilder::RawTowerBuilder(const std::string &name)
  : SubsysReco(name)
{
}

int RawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  std::string geonodename = "CYLINDERCELLGEOM_" + m_Detector;
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geonodename);
  if (!cellgeos)
  {
    std::cout << PHWHERE << " " << geonodename
              << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + geonodename + " node in RawTowerBuilder::CreateNodes");
  }

  // fill the number of layers in the calorimeter
  m_NumLayers = cellgeos->get_NLayers();

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    //exit(1);
  }

  if (Verbosity() >= 1)
  {
    std::cout << "RawTowerBuilder::InitRun :";
    if (m_TowerEnergySrcEnum == kEnergyDeposition)
    {
      std::cout << "save Geant4 energy deposition as the weight of the cells"
                << std::endl;
    }
    else if (m_TowerEnergySrcEnum == kLightYield)
    {
      std::cout << "save light yield as the weight of the cells" << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
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
  std::string cellnodename = "G4CELL_" + m_Detector;
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    std::cout << PHWHERE << " " << cellnodename
              << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // loop over all cells in an event
  PHG4CellContainer::ConstIterator cell_iter;
  PHG4CellContainer::ConstRange cell_range = cells->getCells();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second; ++cell_iter)
  {
    PHG4Cell *cell = cell_iter->second;

    if (Verbosity() > 2)
    {
      std::cout << PHWHERE << " print out the cell:" << std::endl;
      cell->identify();
    }


    //Calculate the cell weight
    float cell_weight = 0;
    if (m_TowerEnergySrcEnum == kEnergyDeposition)
    {
      cell_weight = cell->get_edep();
    }
    else if (m_TowerEnergySrcEnum == kLightYield)
    {
      cell_weight = cell->get_light_yield();
    }


    // add the energy to the corresponding tower
    if (m_UseTowerInfo  != 1)
      {
	RawTower *tower = nullptr;
	short int firstpar;
	short int secondpar;
	if (cell->has_binning(PHG4CellDefs::sizebinning))
	  {
	    firstpar = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());
	    secondpar = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());
	  }
	else if (cell->has_binning(PHG4CellDefs::spacalbinning))
	  {
	    firstpar = PHG4CellDefs::SpacalBinning::get_etabin(cell->get_cellid());
	    secondpar = PHG4CellDefs::SpacalBinning::get_phibin(cell->get_cellid());
	  }
	else
	  {
	    boost::io::ios_flags_saver ifs(std::cout);
	    std::cout << "unknown cell binning, implement 0x" << std::hex << PHG4CellDefs::get_binning(cell->get_cellid()) << std::dec << std::endl;
	    exit(1);
	  }
	tower = m_TowerContainer->getTower(firstpar, secondpar);
	if (!tower)
	  {
	    tower = new RawTowerv1();
	    tower->set_energy(0);
	    m_TowerContainer->AddTower(firstpar, secondpar, tower);
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
	
	if (Verbosity() > 2)
	  {
	    m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
	    tower->identify();
	  }
      }
    if (m_UseTowerInfo > 0)
      {
	TowerInfo *towerinfo;
	unsigned int etabin;
	unsigned int phibin;
	if (cell->has_binning(PHG4CellDefs::spacalbinning))
	  {
	    etabin = PHG4CellDefs::SpacalBinning::get_etabin(cell->get_cellid());
	    phibin = PHG4CellDefs::SpacalBinning::get_phibin(cell->get_cellid());
	  }
	else
	  {
	    boost::io::ios_flags_saver ifs(std::cout);
	    std::cout << "unknown cell binning, implement 0x" << std::hex << PHG4CellDefs::get_binning(cell->get_cellid()) << std::dec << std::endl;
	    exit(1);
	  }
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

  if (m_UseTowerInfo != 1 )
    {
      double towerE = 0;
      if (m_ChkEnergyConservationFlag)
	{
	  double cellE = cells->getTotalEdep();
	  towerE = m_TowerContainer->getTotalEdep();
	  if (fabs(cellE - towerE) / cellE > 1e-5)
	    {
	      std::cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
			<< cellE - towerE << std::endl;
	    }
	}
      if (Verbosity())
	{
	  towerE = m_TowerContainer->getTotalEdep();
	}
      
      m_TowerContainer->compress(m_Emin);
      if (Verbosity())
	{
	  std::cout << "Energy lost by dropping towers with less than " << m_Emin
		    << " GeV energy, lost energy: " << towerE - m_TowerContainer->getTotalEdep()
		    << std::endl;
	  m_TowerContainer->identify();
	  RawTowerContainer::ConstRange begin_end = m_TowerContainer->getTowers();
	  RawTowerContainer::ConstIterator iter;
	  for (iter = begin_end.first; iter != begin_end.second; ++iter)
	    {
	      iter->second->identify();
	    }
	}
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawTowerBuilder::CreateNodes");
  }

  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }

  const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(m_Detector);

  // get the cell geometry and build up the tower geometry object
  std::string geonodename = "CYLINDERCELLGEOM_" + m_Detector;
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geonodename);
  if (!cellgeos)
  {
    std::cout << PHWHERE << " " << geonodename
              << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find " + geonodename + " node in RawTowerBuilder::CreateNodes");
  }
  m_TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!m_RawTowerGeomContainer)
  {
    m_RawTowerGeomContainer = new RawTowerGeomContainer_Cylinderv1(caloid);
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeomContainer, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  // fill the number of layers in the calorimeter
  m_NumLayers = cellgeos->get_NLayers();

  // Create the tower nodes on the tree
  PHG4CylinderCellGeomContainer::ConstIterator miter;
  PHG4CylinderCellGeomContainer::ConstRange begin_end =
      cellgeos->get_begin_end();
  int ifirst = 1;
  int first_layer = -1;
  PHG4CylinderCellGeom *first_cellgeo = nullptr;
  double inner_radius = 0;
  double thickness = 0;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
  {
    PHG4CylinderCellGeom *cellgeo = miter->second;

    if (Verbosity())
    {
      cellgeo->identify();
    }
    thickness += cellgeo->get_thickness();
    if (ifirst)
    {
      first_cellgeo = miter->second;
      m_CellBinning = cellgeo->get_binning();
      m_NumPhiBins = cellgeo->get_phibins();
      m_PhiMin = cellgeo->get_phimin();
      m_PhiStep = cellgeo->get_phistep();
      if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning)
      {
        m_NumEtaBins = cellgeo->get_etabins();
        m_EtaMin = cellgeo->get_etamin();
        m_EtaStep = cellgeo->get_etastep();
      }
      else if (m_CellBinning == PHG4CellDefs::sizebinning)
      {
        m_NumEtaBins = cellgeo->get_zbins();  // bin eta in the same number of z bins
      }
      else if (m_CellBinning == PHG4CellDefs::spacalbinning)
      {
        // use eta definiton for each row of towers
        m_NumEtaBins = cellgeo->get_etabins();
      }
      else
      {
        std::cout << "RawTowerBuilder::CreateNodes::" << Name()
                  << " - Fatal Error - unsupported cell binning method "
                  << m_CellBinning << std::endl;
      }
      inner_radius = cellgeo->get_radius();
      first_layer = cellgeo->get_layer();
      ifirst = 0;
    }
    else
    {
      if (m_CellBinning != cellgeo->get_binning())
      {
        std::cout << "inconsistent binning method from " << m_CellBinning
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_binning() << std::endl;
        exit(1);
      }
      if (inner_radius > cellgeo->get_radius())
      {
        std::cout << "radius of layer " << cellgeo->get_layer() << " is "
                  << cellgeo->get_radius() << " which smaller than radius "
                  << inner_radius << " of first layer in list " << first_layer
                  << std::endl;
      }
      if (m_NumPhiBins != cellgeo->get_phibins())
      {
        std::cout << "mixing number of phibins, fisrt layer: " << m_NumPhiBins
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_phibins() << std::endl;
        exit(1);
      }
      if (m_PhiMin != cellgeo->get_phimin())
      {
        std::cout << "mixing number of phimin, fisrt layer: " << m_PhiMin
                  << " layer " << cellgeo->get_layer() << ": "
                  << cellgeo->get_phimin() << std::endl;
        exit(1);
      }
      if (m_PhiStep != cellgeo->get_phistep())
      {
        std::cout << "mixing phisteps first layer: " << m_PhiStep << " layer "
                  << cellgeo->get_layer() << ": " << cellgeo->get_phistep()
                  << " diff: " << m_PhiStep - cellgeo->get_phistep() << std::endl;
        exit(1);
      }
      if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning)
      {
        if (m_NumEtaBins != cellgeo->get_etabins())
        {
          std::cout << "mixing number of EtaBins , first layer: "
                    << m_NumEtaBins << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_etabins() << std::endl;
          exit(1);
        }
        if (fabs(m_EtaMin - cellgeo->get_etamin()) > 1e-9)
        {
          std::cout << "mixing etamin, fisrt layer: " << m_EtaMin << " layer "
                    << cellgeo->get_layer() << ": " << cellgeo->get_etamin()
                    << " diff: " << m_EtaMin - cellgeo->get_etamin() << std::endl;
          exit(1);
        }
        if (fabs(m_EtaStep - cellgeo->get_etastep()) > 1e-9)
        {
          std::cout << "mixing eta steps first layer: " << m_EtaStep
                    << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_etastep() << " diff: "
                    << m_EtaStep - cellgeo->get_etastep() << std::endl;
          exit(1);
        }
      }

      else if (m_CellBinning == PHG4CellDefs::sizebinning)
      {
        if (m_NumEtaBins != cellgeo->get_zbins())
        {
          std::cout << "mixing number of z bins , first layer: " << m_NumEtaBins
                    << " layer " << cellgeo->get_layer() << ": "
                    << cellgeo->get_zbins() << std::endl;
          exit(1);
        }
      }
    }
  }
  m_RawTowerGeomContainer->set_radius(inner_radius);
  m_RawTowerGeomContainer->set_thickness(thickness);
  m_RawTowerGeomContainer->set_phibins(m_NumPhiBins);
  //  m_RawTowerGeomContainer->set_phistep(m_PhiStep);
  //  m_RawTowerGeomContainer->set_phimin(m_PhiMin);
  m_RawTowerGeomContainer->set_etabins(m_NumEtaBins);

  if (!first_cellgeo)
  {
    std::cout << "RawTowerBuilder::CreateNodes - ERROR - can not find first layer of cells "
              << std::endl;

    exit(1);
  }

  for (int ibin = 0; ibin < first_cellgeo->get_phibins(); ibin++)
  {
    const std::pair<double, double> range = first_cellgeo->get_phibounds(ibin);

    m_RawTowerGeomContainer->set_phibounds(ibin, range);
  }

  if (m_CellBinning == PHG4CellDefs::etaphibinning || m_CellBinning == PHG4CellDefs::etaslatbinning || m_CellBinning == PHG4CellDefs::spacalbinning)
  {
    const double r = inner_radius;

    for (int ibin = 0; ibin < first_cellgeo->get_etabins(); ibin++)
    {
      const std::pair<double, double> range = first_cellgeo->get_etabounds(ibin);

      m_RawTowerGeomContainer->set_etabounds(ibin, range);
    }

    // setup location of all towers
    for (int iphi = 0; iphi < m_RawTowerGeomContainer->get_phibins(); iphi++)
    {
      for (int ieta = 0; ieta < m_RawTowerGeomContainer->get_etabins(); ieta++)
      {
        const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid(caloid, ieta, iphi);

        const double x(r * cos(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double y(r * sin(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double z(r / tan(PHG4Utils::get_theta(m_RawTowerGeomContainer->get_etacenter(ieta))));

        RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
        if (tg)
        {
          if (Verbosity() > 0)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Tower geometry " << key << " already exists" << std::endl;
          }

          if (fabs(tg->get_center_x() - x) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                      << std::endl;

            exit(1);
          }
          if (fabs(tg->get_center_y() - y) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                      << std::endl;
            exit(1);
          }
          if (fabs(tg->get_center_z() - z) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                      << std::endl;
            exit(1);
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << "RawTowerBuilder::CreateNodes - building tower geometry " << key << "" << std::endl;
          }

          tg = new RawTowerGeomv1(key);

          tg->set_center_x(x);
          tg->set_center_y(y);
          tg->set_center_z(z);
          m_RawTowerGeomContainer->add_tower_geometry(tg);
        }
      }
    }
  }
  else if (m_CellBinning == PHG4CellDefs::sizebinning)
  {
    const double r = inner_radius;

    for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
    {
      const std::pair<double, double> z_range = first_cellgeo->get_zbounds(ibin);
      //          const double r = first_cellgeo->get_radius();
      const double eta1 = -log(tan(atan2(r, z_range.first) / 2.));
      const double eta2 = -log(tan(atan2(r, z_range.second) / 2.));
      m_RawTowerGeomContainer->set_etabounds(ibin, std::make_pair(eta1, eta2));
    }

    // setup location of all towers
    for (int iphi = 0; iphi < m_RawTowerGeomContainer->get_phibins(); iphi++)
    {
      for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(caloid, ibin, iphi);

        const double x(r * cos(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double y(r * sin(m_RawTowerGeomContainer->get_phicenter(iphi)));
        const double z(first_cellgeo->get_zcenter(ibin));

        RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
        if (tg)
        {
          if (Verbosity() > 0)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Tower geometry " << key << " already exists" << std::endl;
          }

          if (fabs(tg->get_center_x() - x) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                      << std::endl;

            exit(1);
          }
          if (fabs(tg->get_center_y() - y) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                      << std::endl;
            exit(1);
          }
          if (fabs(tg->get_center_z() - z) > 1e-4)
          {
            std::cout << "RawTowerBuilder::CreateNodes - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                      << std::endl;
            exit(1);
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            std::cout << "RawTowerBuilder::CreateNodes - building tower geometry " << key << "" << std::endl;
          }

          tg = new RawTowerGeomv1(key);

          tg->set_center_x(x);
          tg->set_center_y(y);
          tg->set_center_z(z);
          m_RawTowerGeomContainer->add_tower_geometry(tg);
        }
      }
    }
  }
  else
  {
    std::cout << "RawTowerBuilder::CreateNodes - ERROR - unsupported cell geometry "
              << m_CellBinning << std::endl;
    exit(1);
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeomContainer->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerBuilder::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  
  if (m_UseTowerInfo != 1)
    {
      // Create the tower nodes on the tree
      if (m_SimTowerNodePrefix.empty())
	{
	  // no prefix, consistent with older convention
	  m_TowerNodeName = "TOWER_" + m_Detector;
	}
      else
	{
	  m_TowerNodeName = "TOWER_" + m_SimTowerNodePrefix + "_" + m_Detector;
	}
      m_TowerContainer = findNode::getClass<RawTowerContainer>(DetNode, m_TowerNodeName.c_str());
      if (m_TowerContainer == nullptr)
	{
	  m_TowerContainer = new RawTowerContainer(caloid);
	  
	  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_TowerContainer, m_TowerNodeName, "PHObject");
	  DetNode->addNode(towerNode);
	}
    }



  if (m_UseTowerInfo > 0 )
    {
     if (m_SimTowerNodePrefix.empty())
	{
	  m_TowerInfoNodeName = "TOWERINFO_" + m_Detector;
	}
      else
	{
	  m_TowerInfoNodeName = "TOWERINFO_" + m_SimTowerNodePrefix + "_" + m_Detector;
	}
     TowerInfoContainer* m_TowerInfoContainer = findNode::getClass<TowerInfoContainer>(DetNode,m_TowerInfoNodeName);
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
