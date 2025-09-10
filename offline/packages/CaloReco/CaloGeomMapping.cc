#include "CaloGeomMapping.h"

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
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerGeomv5.h>

#include <cmath>      // for fabs, atan, cos
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, endl
#include <stdexcept>  // for runtime_error
#include <utility>    // for pair

//____________________________________________________________________________..
CaloGeomMapping::CaloGeomMapping(const std::string &name)
  : SubsysReco(name)
  , m_Detector("CEMC")
{
}

//____________________________________________________________________________..
int CaloGeomMapping::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "CaloGeomMapping::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  }
  /* std::cout << "Printing node tree before new node creation:" << std::endl; */
  /* topNode->print(); */
  try
  {
    CreateGeomNode(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(1);
  }
  /* std::cout << "Printing node tree after new node creation:" << std::endl; */
  /* topNode->print(); */
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloGeomMapping::CreateGeomNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in CaloGeomMapping::CreateGeomNode");
  }

  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", m_Detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(RunDetNode);
  }

  m_caloid = RawTowerDefs::convert_name_to_caloid(m_Detector);
  m_TowerGeomNodeName = "TOWERGEOM_" + m_Detector;

  if (m_UseDetailedGeometry && m_Detector == "CEMC")
  {
    // Rename the TOWERGEOM node to avoid confusion with the former geometry
    m_TowerGeomNodeName = m_TowerGeomNodeName + "_DETAILED";
  }

  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, m_TowerGeomNodeName);
  if (!m_RawTowerGeomContainer)
  {
    m_RawTowerGeomContainer = new RawTowerGeomContainer_Cylinderv1(m_caloid);
    // add it to the node tree
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeomContainer, m_TowerGeomNodeName, "PHObject");
    RunDetNode->addNode(newNode);
  }
  
  // Get the geometry mapping file from the Conditions Database
  std::string inName=CDBInterface::instance()->getUrl("CALO_TOWER_GEOMETRY");
  CDBTTree * cdbttree = new CDBTTree(inName);
  cdbttree->LoadCalibrations();

  std::string parName;
  std::string parBase;
  // Set the radius, thickness, number of eta and phi bins
  if (m_Detector == "CEMC")
  {
    parBase = "cemc";
    m_RawTowerGeomContainer->set_radius(93.5);
    m_RawTowerGeomContainer->set_thickness(20.4997);
    m_RawTowerGeomContainer->set_phibins(256);
    m_RawTowerGeomContainer->set_etabins(96);
  }
  if (m_Detector == "HCALIN")
  {
    parBase = "hcalin";
    m_RawTowerGeomContainer->set_radius(115);
    m_RawTowerGeomContainer->set_thickness(25.005);
    m_RawTowerGeomContainer->set_phibins(64);
    m_RawTowerGeomContainer->set_etabins(24);
  }
  if (m_Detector == "HCALOUT")
  {
    parBase = "hcalout";
    m_RawTowerGeomContainer->set_radius(177.423);
    m_RawTowerGeomContainer->set_thickness(96.894);
    m_RawTowerGeomContainer->set_phibins(64);
    m_RawTowerGeomContainer->set_etabins(24);
  }

  // Set the eta and phi bounds of each bin
  for (int ibin = 0; ibin < m_RawTowerGeomContainer->get_etabins(); ibin++)
  {
    parName = parBase + "_eta_";
    double first;
    double second;
    first = cdbttree->GetDoubleValue(ibin, parName + "first");
    second = cdbttree->GetDoubleValue(ibin, parName + "second");
    const std::pair<double, double> range(first, second);
    m_RawTowerGeomContainer->set_etabounds(ibin, range);
  }
  for (int ibin = 0; ibin < m_RawTowerGeomContainer->get_phibins(); ibin++)
  {
    parName = parBase + "_phi_";
    double first;
    double second;
    first = cdbttree->GetDoubleValue(ibin, parName + "first");
    second = cdbttree->GetDoubleValue(ibin, parName + "second");
    const std::pair<double, double> range(first, second);
    m_RawTowerGeomContainer->set_phibounds(ibin, range);
  }


  // Build the RawTowerGeom objects
  if (m_UseDetailedGeometry == true && m_Detector == "CEMC")
  {
    // Detailed geometry using the 8 block vertices as defined in the GEANT4 simulation
    BuildDetailedGeometry();
  }
  else 
  {
    // Approximate geometry
    // Deduces the tower position from the eta/phi bins
    // and project them at a constant radius.
    if (m_UseDetailedGeometry == true)
    {
      std::cout << "CaloGeomMapping::CreateNodes - Detailed geometry is not yet supported for " << m_Detector << ". The former approximate geometry is used instead." << std::endl;
    }
    BuildFormerGeometry();
  }
}

void CaloGeomMapping::BuildFormerGeometry()
{
  // Populate container with RawTowerGeom objects
  for (int ieta = 0; ieta < m_RawTowerGeomContainer->get_etabins(); ieta++)
  {
    for (int iphi = 0; iphi < m_RawTowerGeomContainer->get_phibins(); iphi++)
    {
      // build tower geom here
      const RawTowerDefs::keytype key =
          RawTowerDefs::encode_towerid(m_caloid, ieta, iphi);

      double r = m_RawTowerGeomContainer->get_radius();
      const double x(r * cos(m_RawTowerGeomContainer->get_phicenter(iphi)));
      const double y(r * sin(m_RawTowerGeomContainer->get_phicenter(iphi)));
      const double z(r / tan(2 * atan(exp(-1 * m_RawTowerGeomContainer->get_etacenter(ieta)))));
      RawTowerGeom *tg = m_RawTowerGeomContainer->get_tower_geometry(key);
      if (tg)
      {
        if (Verbosity() > 0)
        {
          std::cout << "CaloGeomMapping::BuildFormerGeometry - Tower geometry " << key << " already exists" << std::endl;
        }

        if (fabs(tg->get_center_x() - x) > 1e-4)
        {
          std::cout << "CaloGeomMapping::BuildFormerGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                    << std::endl;

          exit(1);
        }
        if (fabs(tg->get_center_y() - y) > 1e-4)
        {
          std::cout << "CaloGeomMapping::BuildFormerGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                    << std::endl;
          exit(1);
        }
        if (fabs(tg->get_center_z() - z) > 1e-4)
        {
          std::cout << "CaloGeomMapping::BuildFormerGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                    << std::endl;
          exit(1);
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          std::cout << "CaloGeomMapping::BuildFormerGeometry - building tower geometry " << key << "" << std::endl;
        }

        tg = new RawTowerGeomv1(key);

        tg->set_center_x(x);
        tg->set_center_y(y);
        tg->set_center_z(z);
        m_RawTowerGeomContainer->add_tower_geometry(tg);
      }
    }
  }  // end loop over eta, phi bins
}

void CaloGeomMapping::BuildDetailedGeometry()
{
  // This method is only implemented for the electromagnetic calorimeter

  // Get the geometry mapping file from the Conditions Database
  std::string inName = CDBInterface::instance()->getUrl("CALO_TOWER_GEOMETRY_EMCAL_DETAILED");
  CDBTTree *cdbttree = new CDBTTree(inName);
  cdbttree->LoadCalibrations();
  //std::cout << "printing CDB TTree" << std::endl;
  //cdbttree->Print();

  // Populate container with RawTowerGeom objects
  const int nSectors = 32;
  const int nUniqueTowers = 192;
  for (int iSector = 0; iSector < nSectors; iSector++)
  {
    for (int iUniqueTower = 0; iUniqueTower < nUniqueTowers; iUniqueTower++)
    {
      // Vertices' coordinates of the corresponding unique block
      // separated into 4 towers (2 bins in eta x 2 bins in phi)
      
      /*   R                     
           |                    7_______14      14_______6
           |                   /|      /|       /|      /|
           |_______ eta     15/_|_____/18    18/_|_____/13
          /                   | |     | |      | |     | |
         /                    |3|_____|_|10    10|_____|_|2
        /                     | /     | /      | /     | /
       phi                  11|/______|/17   17|/______|/9

                          15_______18      18_______13
                          /|      /|       /|      /|
                        8/_|_____/16    16/_|_____/5| 
                         | |     | |      | |     | |
                         11|_____|_|17    17|_____|_|9 
                         | /     | /      | /     | /
                        4|/______|/12   12|/______|/1
      */


      std::vector<double> vertices(18*3);
      // Block vertices (1 -> 8)
      for (int i = 0; i < 8; i++) 
      {
        vertices[i*3+0] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_x") / 10; // Convert mm to cm
        vertices[i*3+1] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_y") / 10;
        vertices[i*3+2] = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1, "vtx_" + std::to_string(i) + "_z") / 10;
      }

      // Get the rotation of the tower:
      double rot_x = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1,"rot_x");
      double rot_y = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1,"rot_y");
      double rot_z = cdbttree->GetDoubleValue(iUniqueTower * nSectors + iSector + 1,"rot_z");

      // Divide each block into 4 equal towers.
      for (int j = 0; j < 3; j++) 
      {
        for (int i = 0; i < 4; i++)
        {
          vertices[(8 + i) * 3 + j]=(vertices[i * 3 + j] + vertices[((i + 1) % 4) * 3 + j]) / 2; // Middle points corresponding to tower vertices (inner side)
          vertices[(8 + i + 4) * 3 + j]=(vertices[(i + 4) * 3 + j] + vertices[((i + 1) % 4 + 4) * 3 + j]) / 2; // Middle points (outer side)
        }
        vertices[(16)*3+j] = (vertices[0 * 3 + j] + vertices[1 * 3 + j] + vertices[2 * 3 + j] + vertices[3 * 3 + j]) / 4;
        vertices[(17)*3+j] = (vertices[4 * 3 + j] + vertices[5 * 3 + j] + vertices[6 * 3 + j] + vertices[7 * 3 + j]) / 4;
      }

      int ietaBlock = iUniqueTower / 4;
      int iphiBlock = iSector * 4 + iUniqueTower % 4; 

      // Vertices' indices for each of the 4 individual towers.
      double sub_tower[2][2][8];
      double sub_tower_01[8] = {12, 17, 11, 4, 16, 18, 15, 8};
      double sub_tower_00[8] = {17, 10, 3, 11, 18, 14, 7, 15};
      double sub_tower_11[8] = {1, 9, 17, 12, 5, 13, 18, 16};
      double sub_tower_10[8] = {9, 2, 10, 17, 13, 6, 14, 18};
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

          std::vector<double> towerVertices(8*3);
          for (int ivtx = 0; ivtx < 8; ivtx++)
          {
            for (int icoord = 0; icoord < 3; icoord++) 
            {
              towerVertices[ivtx*3+icoord] = vertices[(sub_tower[iTower][jTower][ivtx]-1)*3+icoord];
            }
          }
          
          // build tower geom here
          const RawTowerDefs::keytype key =
            RawTowerDefs::encode_towerid(m_caloid, ieta, iphi);
          
      
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
              std::cout << "CaloGeomMapping::BuildDetailedGeometry - Tower geometry " << key << " already exists" << std::endl;
            }

            if (fabs(tg->get_center_x() - x) > 1e-4)
            {
              std::cout << "CaloGeomMapping::BuildDetailedGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing x = " << tg->get_center_x() << " and expected x = " << x
                        << std::endl;
              
              exit(1);
            }
            if (fabs(tg->get_center_y() - y) > 1e-4)
            {
              std::cout << "CaloGeomMapping::BuildDetailedGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing y = " << tg->get_center_y() << " and expected y = " << y
                        << std::endl;
              exit(1);
            }
            if (fabs(tg->get_center_z() - z) > 1e-4)
            {
              std::cout << "CaloGeomMapping::BuildDetailedGeometry - Fatal Error - duplicated Tower geometry " << key << " with existing z= " << tg->get_center_z() << " and expected z = " << z
                        << std::endl;
              exit(1);
            }
          }
          else
          {
            if (Verbosity() > 0)
            {
              std::cout << "CaloGeomMapping::BuildDetailedGeometry - building tower geometry " << key << "" << std::endl;
            }
            
            tg = tg0;
            m_RawTowerGeomContainer->add_tower_geometry(tg);
          }
        }
      }
    }
  }  // end loop over eta, phi bins
}

