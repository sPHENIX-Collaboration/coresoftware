#include "CaloGeomMappingv2.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <calobase/RawTowerDefs.h>           // for encode_towerid
#include <calobase/RawTowerGeom.h>           // for RawTowerGeom
#include <calobase/RawTowerGeomContainer.h>  // for RawTowerGeomC...
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv5.h>

#include <cmath>      // for fabs, atan, cos
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, endl
#include <stdexcept>  // for runtime_error
#include <utility>    // for pair

//____________________________________________________________________________..
CaloGeomMappingv2::CaloGeomMappingv2(const std::string &name)
  : SubsysReco(name)
  , m_Detector("CEMC")
  , m_RawTowerGeomContainer(nullptr)
{
  std::cout << "CaloGeomMappingv2::CaloGeomMappingv2(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloGeomMappingv2::~CaloGeomMappingv2()
{
  std::cout << "CaloGeomMappingv2::~CaloGeomMappingv2() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloGeomMappingv2::Init(PHCompositeNode *topNode)
{
  std::cout << "CaloGeomMappingv2::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  try
  {
    CreateGeomNode(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloGeomMappingv2::CreateGeomNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in CaloGeomMappingv2::CreateGeomNode");
  }

  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }

  const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(m_Detector);
  if (m_TowerGeomNodeName.empty())
  {
    m_TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  }
  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!m_RawTowerGeomContainer)
  {
    m_RawTowerGeomContainer = new RawTowerGeomContainer_Cylinderv1(caloid);
    // add it to the node tree
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeomContainer, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }

  // Get the geometry mapping file from the Conditions Database
  // For the moment, the default geometry file is not yet in CDB

  // All towers' geome
  std::string inName = "/sphenix/user/virgilemahaut/geometry/macros/calo_geom_mapping_exact.root";
  CDBTTree *cdbttree = new CDBTTree(inName);

  std::string parBase;
  // Set the radius, thickness, number of eta and phi bins
  if (m_Detector == "CEMC")
  {
    parBase = "cemc";
    m_RawTowerGeomContainer->set_radius(93.5);
    m_RawTowerGeomContainer->set_thickness(20.4997);
    m_RawTowerGeomContainer->set_phibins(256);
    m_RawTowerGeomContainer->set_etabins(96);
    //  m_RawTowerGeomContainer->set_phistep(m_PhiStep);
    //  m_RawTowerGeomContainer->set_phimin(m_PhiMin);
  }
  if (m_Detector == "HCALIN")
  {
    parBase = "hcalin";
    m_RawTowerGeomContainer->set_radius(115);
    m_RawTowerGeomContainer->set_thickness(25.005);
    m_RawTowerGeomContainer->set_phibins(64);
    m_RawTowerGeomContainer->set_etabins(24);
    //  m_RawTowerGeomContainer->set_phistep(m_PhiStep);
    //  m_RawTowerGeomContainer->set_phimin(m_PhiMin);
  }
  if (m_Detector == "HCALOUT")
  {
    parBase = "hcalout";
    m_RawTowerGeomContainer->set_radius(177.423);
    m_RawTowerGeomContainer->set_thickness(96.894);
    m_RawTowerGeomContainer->set_phibins(64);
    m_RawTowerGeomContainer->set_etabins(24);
    //  m_RawTowerGeomContainer->set_phistep(m_PhiStep);
    //  m_RawTowerGeomContainer->set_phimin(m_PhiMin);
  }

  // Populate container with RawTowerGeom objects
  // The calorimeter is built out of 32 sectors with 192 unique towers (here, a "unique tower" actually corresponds to a block
  const int nSectors = 32;
  const int nUniqueTowers = 192; 
  const int nDimensions = 3;
  // Each block is made of 4 neighbouring towers, for a total of 18 separate tower vertices.
  const int nTowerVerticesPerBlock = 18;
  const int nVerticesPerTower = 8;
  for (int iSector = 0; iSector < nSectors; iSector++)
  {
    for (int iUniqueTower = 0; iUniqueTower < nUniqueTowers; iUniqueTower++)
    {
      // Vertices' coordinates of the corresponding unique block.
      std::vector<double> vertices((unsigned long) nTowerVerticesPerBlock * nDimensions);
      for (int i = 0; i < 8; i++)
      {
        vertices[i * nDimensions + 0] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_x") / 10;  // Convert mm to cm
        vertices[i * nDimensions + 1] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_y") / 10;
        vertices[i * nDimensions + 2] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_z") / 10;
      }

      // Get the rotation vector of the block:
      double rot_x = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "rot_x");
      double rot_y = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "rot_y");
      double rot_z = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "rot_z");

      // Divide each block into 4 equal towers.
      for (int j = 0; j < 3; j++)
      {
        for (int i = 0; i < 4; i++)
        {
          vertices[(8 + i) * 3 + j] = (vertices[i * 3 + j] + vertices[((i + 1) % 4) * 3 + j]) / 2;
          vertices[(8 + i + 4) * 3 + j] = (vertices[(i + 4) * 3 + j] + vertices[((i + 1) % 4 + 4) * 3 + j]) / 2;
        }
        vertices[(16) * 3 + j] = (vertices[0 * 3 + j] + vertices[1 * 3 + j] + vertices[2 * 3 + j] + vertices[3 * 3 + j]) / 4;
        vertices[(17) * 3 + j] = (vertices[4 * 3 + j] + vertices[5 * 3 + j] + vertices[6 * 3 + j] + vertices[7 * 3 + j]) / 4;
      }

      int ietaBlock = iUniqueTower / 4;
      int iphiBlock = iSector * 4 + iUniqueTower % 4;

      // Vertices' indices for each of the 4 individual towers.
      double sub_tower[2][2][8];
      double sub_tower_01[8] = {12, 17, 11, 4, 16, 18, 15, 8};  // low eta high phi
      double sub_tower_00[8] = {17, 10, 3, 11, 18, 14, 7, 15};  // low eta low phi
      double sub_tower_11[8] = {1, 9, 17, 12, 5, 13, 18, 16};   // high eta high phi
      double sub_tower_10[8] = {9, 2, 10, 17, 13, 6, 14, 18};   // high eta low phi
      for (int ivtx = 0; ivtx < 8; ivtx++)
      {
        sub_tower[0][0][ivtx] = sub_tower_00[ivtx];
        sub_tower[0][1][ivtx] = sub_tower_01[ivtx];
        sub_tower[1][0][ivtx] = sub_tower_10[ivtx];
        sub_tower[1][1][ivtx] = sub_tower_11[ivtx];
      }

      for (int iTower = 0; iTower < 2; iTower++)
      {
        for (int jTower = 0; jTower < 2; jTower++)
        {
          int ieta = ietaBlock * 2 + iTower;
          int iphi = iphiBlock * 2 + jTower;

          std::vector<double> towerVertices((unsigned long) nVerticesPerTower * nDimensions); // 8 x 3
          for (int ivtx = 0; ivtx < nVerticesPerTower; ivtx++)
          {
            for (int icoord = 0; icoord < nDimensions; icoord++)
            {
              towerVertices[ivtx * nDimensions + icoord] = vertices[(sub_tower[iTower][jTower][ivtx] - 1) * nDimensions + icoord];
            }
          }

          // build tower geom here
          const RawTowerDefs::keytype key =
              RawTowerDefs::encode_towerid(caloid, ieta, iphi);

          RawTowerGeomv5 *tg0 = new RawTowerGeomv5(key);

          tg0->set_vertices(towerVertices);
          tg0->set_rotx(rot_x);
          tg0->set_roty(rot_y);
          tg0->set_rotz(rot_z);

          const double x(tg0->get_center_x());
          const double y(tg0->get_center_y());
          const double z(tg0->get_center_z());

          RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
          if (tg)
          {
            delete tg0;
            if (Verbosity() > 0)
            {
              std::cout << "CaloGeomMappingv2::CreateGeomNode - Tower geometry " << key << " already exists" << std::endl;
            }

            if (fabs(tg->get_center_x() - x) > 1e-4)
            {
              std::cout << "CaloGeomMappingv2::CreateGeomNode - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                        << std::endl;

              exit(1);
            }
            if (fabs(tg->get_center_y() - y) > 1e-4)
            {
              std::cout << "CaloGeomMappingv2::CreateGeomNode - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                        << std::endl;
              exit(1);
            }
            if (fabs(tg->get_center_z() - z) > 1e-4)
            {
              std::cout << "CaloGeomMappingv2::CreateGeomNode - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                        << std::endl;
              exit(1);
            }
          }
          else
          {
            if (Verbosity() > 0)
            {
              std::cout << "CaloGeomMappingv2::CreateGeomNode - building tower geometry " << key << "" << std::endl;
            }

            tg = tg0;
            m_RawTowerGeomContainer->add_tower_geometry(tg);
          }
        }
      }
    }
  }  // end loop over eta, phi bins
}  // end of building RawTowerGeomContainer

void CaloGeomMappingv2::set_detector_name(const std::string &name)
{
  m_Detector = name;
}

void CaloGeomMappingv2::setTowerGeomNodeName(const std::string &name)
{
  m_TowerGeomNodeName = name;
}

std::string CaloGeomMappingv2::get_detector_name()
{
  return m_Detector;
}
