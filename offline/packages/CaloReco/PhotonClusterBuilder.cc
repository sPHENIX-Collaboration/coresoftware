#include "PhotonClusterBuilder.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawCluster.h>
#include <calobase/PhotonClusterContainer.h>
#include <calobase/PhotonClusterv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeomContainer.h>


// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>
#include <phool/PHNodeIterator.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TMVA/RBDT.hxx>

#include <iostream>
#include <stdexcept>
#include <cmath>

namespace
{
  // Helper function to shift tower indices for wrapping in phi
  void shift_tower_index(int& ieta, int& iphi, int etadiv, int phidiv)
  {
    while (iphi < 0) iphi += phidiv;
    while (iphi >= phidiv) iphi -= phidiv;
    if (ieta < 0 || ieta >= etadiv)
    {
      ieta = -1; // invalid
    }
  }
}


PhotonClusterBuilder::PhotonClusterBuilder(const std::string& name)
  : SubsysReco(name)
{
}

int PhotonClusterBuilder::InitRun(PHCompositeNode* topNode)
{
  // BDT
  m_bdt = std::make_unique<TMVA::Experimental::RBDT>("myBDT", m_bdt_model_file);

  // locate input raw cluster container
  m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
  if (!m_rawclusters)
  {
    std::cerr << Name() << ": could not find RawClusterContainer node '" << m_input_cluster_node << "'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_emc_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (!m_emc_tower_container)
  {
      std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_CEMC'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!m_geomEM)
  {
      std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_CEMC'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ihcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  if (!m_ihcal_tower_container)
  {
      std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALIN'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!m_geomIH)
  {
      std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALIN'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ohcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if (!m_ohcal_tower_container)
  {
      std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALOUT'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!m_geomOH)
  {
      std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALOUT'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }


  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    throw std::runtime_error("PhotonClusterBuilder: DST node not found");
  }
  m_photon_container = findNode::getClass<PhotonClusterContainer>(dstNode, m_output_photon_node);
  if (!m_photon_container)
  {
    m_photon_container = new PhotonClusterContainer();
    auto photonNode = new PHIODataNode<PHObject>(m_photon_container, m_output_photon_node, "PHObject");
    dstNode->addNode(photonNode);
  }
}

int PhotonClusterBuilder::process_event(PHCompositeNode* topNode)
{
  if (!m_rawclusters)
  {
    m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
    if (!m_rawclusters)
    {
      std::cerr << Name() << ": missing RawClusterContainer '" << m_input_cluster_node << "'" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  //init with NaN
  m_vertex = std::numeric_limits<float>::quiet_NaN();
    //assume we need vertex for photon shower shape for now
    //in the future we need to change the vertex to MBD tracking combined
    MbdVertexMap *vertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");

    if (!vertexmap)
    {
      std::cout << "GlobalVertexMap node is missing" << std::endl;
    }
    if (vertexmap && !vertexmap->empty())
    {
      MbdVertex *vtx = vertexmap->begin()->second;
      if (vtx)
      {
        m_vertex = vtx->get_z();


        if (m_vertex != m_vertex)
          return Fun4AllReturnCodes::EVENT_OK;

      }
      else
      {
        return Fun4AllReturnCodes::EVENT_OK;
      }
    }
  // iterate over clusters via map to have access to keys if needed
  const auto & rcmap = m_rawclusters->getClustersMap();
  for (const auto & kv : rcmap)
  {
    RawCluster* rc = kv.second;
    if (!rc) 
    { 
        continue; 
    }
    
    CLHEP::Hep3Vector vertex_vec(0, 0, m_vertex);

    float eta = RawClusterUtility::GetPseudorapidity(*rc,  vertex_vec);
    float phi = RawClusterUtility::GetAzimuthAngle(*rc,  vertex_vec);
    float E = rc->get_energy();
    float ET = E / cosh(eta);
    if (ET < m_min_cluster_et) 
    { 
        continue; 
    }
    RawClusterv1* rcv1 = dynamic_cast<RawClusterv1*>(rc);
    if (!rcv1) 
    {
        std::cerr << "PhotonClusterBuilder: unsupported RawCluster type" << std::endl;
        continue;
    }
    PhotonClusterv1* photon = new PhotonClusterv1(*rcv1);

    calculate_shower_shapes(rc, photon, eta, phi);
    calculate_bdt_score(photon);

    m_photon_container->AddCluster(photon);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::calculate_bdt_score(PhotonClusterv1* photon)
{
  if (!m_bdt)
  {
    return;
  }

  std::vector<float> x;
  for (const auto& feature : m_bdt_feature_list)
  {
    if (feature == "e11_over_e33")
    {
      float e11 = photon->get_shower_shape_parameter("e11");
      float e33 = photon->get_shower_shape_parameter("e33");
      x.push_back((e33 > 0) ? e11 / e33 : 0);
    }
    else if (feature == "vertex_z")
    {
      x.push_back(m_vertex);
    }
    else
    {
      x.push_back(photon->get_shower_shape_parameter(feature));
    }
  }

  float bdt_score = -1; // default value

  bdt_score = m_bdt->Compute(x)[0];
  photon->set_shower_shape_parameter("bdt_score", bdt_score);
}


void PhotonClusterBuilder::calculate_shower_shapes(RawCluster* rc, PhotonClusterv1* photon, float cluster_eta, float cluster_phi)
{
    std::vector<float> showershape = rc->get_shower_shapes(m_shape_min_tower_E);
    if (showershape.empty())
    {
        return;
    }

    std::pair<int, int> leadtowerindex = rc->get_lead_tower();
    int lead_ieta = leadtowerindex.first;
    int lead_iphi = leadtowerindex.second;

    float avg_eta = showershape[4] + 0.5F;
    float avg_phi = showershape[5] + 0.5F;

    int maxieta = std::floor(avg_eta);
    int maxiphi = std::floor(avg_phi);

    if (maxieta < 3 || maxieta > 92) return;

    // for detamax, dphimax, nsaturated
    int detamax = 0;
    int dphimax = 0;
    int nsaturated = 0;
    const RawCluster::TowerMap tower_map = rc->get_towermap();
    std::set<unsigned int> towers_in_cluster;
    for (auto tower_iter : tower_map)
    {
        towers_in_cluster.insert(tower_iter.first);
        RawTowerDefs::keytype tower_key = tower_iter.first;
        int ieta = RawTowerDefs::decode_index1(tower_key);
        int iphi = RawTowerDefs::decode_index2(tower_key);
        
        unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
        TowerInfo *towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
        if (towerinfo && towerinfo->get_isSaturated())
        {
            nsaturated++;
        }

        int totalphibins = 256;
        auto dphiwrap = [totalphibins](int towerphi, int maxiphi_arg)
        {
          int idphi = towerphi - maxiphi_arg;
          if (idphi > totalphibins / 2) idphi -= totalphibins;
          if (idphi < -totalphibins / 2) idphi += totalphibins;
          return idphi;
        };

        int deta = ieta - lead_ieta;
        int dphi_val = dphiwrap(iphi, lead_iphi);

        if (std::abs(deta) > detamax) detamax = std::abs(deta);
        if (std::abs(dphi_val) > dphimax) dphimax = std::abs(dphi_val);
    }


    float E77[7][7] = {{0.0f}};
    int E77_ownership[7][7] = {{0}};

    for (int ieta = maxieta - 3; ieta < maxieta + 4; ieta++)
    {
        for (int iphi = maxiphi - 3; iphi < maxiphi + 4; iphi++)
        {
            int temp_ieta = ieta;
            int temp_iphi = iphi;
            shift_tower_index(temp_ieta, temp_iphi, 96, 256);
            if (temp_ieta < 0) continue;

            unsigned int towerinfokey = TowerInfoDefs::encode_emcal(temp_ieta, temp_iphi);
            
            if (towers_in_cluster.count(towerinfokey))
            {
                E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 1;
            }
            
            TowerInfo *towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
            if (towerinfo && towerinfo->get_isGood())
            {
                float energy = towerinfo->get_energy();
                if (energy > m_shape_min_tower_E)
                {
                    E77[ieta - maxieta + 3][iphi - maxiphi + 3] = energy;
                }
            }
        }
    }

    float e11 = E77[3][3];
    float e33 = 0, e55 = 0, e77 = 0;
    float e13 = 0, e15 = 0, e17 = 0;
    float e31 = 0, e51 = 0, e71 = 0;
    float e35 = 0, e37 = 0, e53 = 0, e73 = 0, e57 = 0, e75 = 0;
    float weta = 0, wphi = 0, weta_cog = 0, wphi_cog = 0, weta_cogx = 0, wphi_cogx = 0;
    float Eetaphi = 0;
    float shift_eta = avg_eta - std::floor(avg_eta) - 0.5;
    float shift_phi = avg_phi - std::floor(avg_phi) - 0.5;
    float cog_eta = 3 + shift_eta;
    float cog_phi = 3 + shift_phi;
    
    float w32 = 0, e32 = 0, w52 = 0, e52 = 0, w72 = 0, e72 = 0;
    int signphi = (avg_phi - std::floor(avg_phi)) > 0.5 ? 1 : -1;


    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            int di = std::abs(i - 3);
            int dj = std::abs(j - 3);
            float di_float = i - cog_eta;
            float dj_float = j - cog_phi;

            if (E77_ownership[i][j] == 1)
            {
                weta += E77[i][j] * di * di;
                wphi += E77[i][j] * dj * dj;
                weta_cog += E77[i][j] * di_float * di_float;
                wphi_cog += E77[i][j] * dj_float * dj_float;
                Eetaphi += E77[i][j];
                if (i != 3 || j != 3)
                {
                    weta_cogx += E77[i][j] * di_float * di_float;
                    wphi_cogx += E77[i][j] * dj_float * dj_float;
                }
            }

            e77 += E77[i][j];
            if (di <= 1 && (dj == 0 || j == (3 + signphi))) { w32 += E77[i][j] * (i - 3) * (i - 3); e32 += E77[i][j]; }
            if (di <= 2 && (dj == 0 || j == (3 + signphi))) { w52 += E77[i][j] * (i - 3) * (i - 3); e52 += E77[i][j]; }
            if (di <= 3 && (dj == 0 || j == (3 + signphi))) { w72 += E77[i][j] * (i - 3) * (i - 3); e72 += E77[i][j]; }

            if (di <= 0 && dj <= 1) e13 += E77[i][j];
            if (di <= 0 && dj <= 2) e15 += E77[i][j];
            if (di <= 0 && dj <= 3) e17 += E77[i][j];
            if (di <= 1 && dj <= 0) e31 += E77[i][j];
            if (di <= 2 && dj <= 0) e51 += E77[i][j];
            if (di <= 3 && dj <= 0) e71 += E77[i][j];
            if (di <= 1 && dj <= 1) e33 += E77[i][j];
            if (di <= 1 && dj <= 2) e35 += E77[i][j];
            if (di <= 1 && dj <= 3) e37 += E77[i][j];
            if (di <= 2 && dj <= 1) e53 += E77[i][j];
            if (di <= 3 && dj <= 1) e73 += E77[i][j];
            if (di <= 2 && dj <= 2) e55 += E77[i][j];
            if (di <= 2 && dj <= 3) e57 += E77[i][j];
            if (di <= 3 && dj <= 2) e75 += E77[i][j];
        }
    }

    if (Eetaphi > 0)
    {
        weta /= Eetaphi;
        wphi /= Eetaphi;
        weta_cog /= Eetaphi;
        wphi_cog /= Eetaphi;
        weta_cogx /= Eetaphi;
        wphi_cogx /= Eetaphi;
    }
    if (e32 > 0) w32 /= e32;
    if (e52 > 0) w52 /= e52;
    if (e72 > 0) w72 /= e72;


    photon->set_shower_shape_parameter("et1", showershape[0]);
    photon->set_shower_shape_parameter("et2", showershape[1]);
    photon->set_shower_shape_parameter("et3", showershape[2]);
    photon->set_shower_shape_parameter("et4", showershape[3]);
    photon->set_shower_shape_parameter("e11", e11);
    photon->set_shower_shape_parameter("e33", e33);
    photon->set_shower_shape_parameter("e55", e55);
    photon->set_shower_shape_parameter("e77", e77);
    photon->set_shower_shape_parameter("e13", e13);
    photon->set_shower_shape_parameter("e15", e15);
    photon->set_shower_shape_parameter("e17", e17);
    photon->set_shower_shape_parameter("e31", e31);
    photon->set_shower_shape_parameter("e51", e51);
    photon->set_shower_shape_parameter("e71", e71);
    photon->set_shower_shape_parameter("e35", e35);
    photon->set_shower_shape_parameter("e37", e37);
    photon->set_shower_shape_parameter("e53", e53);
    photon->set_shower_shape_parameter("e73", e73);
    photon->set_shower_shape_parameter("e57", e57);
    photon->set_shower_shape_parameter("e75", e75);
    photon->set_shower_shape_parameter("weta", weta);
    photon->set_shower_shape_parameter("wphi", wphi);
    photon->set_shower_shape_parameter("weta_cog", weta_cog);
    photon->set_shower_shape_parameter("wphi_cog", wphi_cog);
    photon->set_shower_shape_parameter("weta_cogx", weta_cogx);
    photon->set_shower_shape_parameter("wphi_cogx", wphi_cogx);
    photon->set_shower_shape_parameter("detamax", detamax);
    photon->set_shower_shape_parameter("dphimax", dphimax);
    photon->set_shower_shape_parameter("nsaturated", nsaturated);
    photon->set_shower_shape_parameter("e32", e32);
    photon->set_shower_shape_parameter("e52", e52);
    photon->set_shower_shape_parameter("e72", e72);
    photon->set_shower_shape_parameter("w32", w32);
    photon->set_shower_shape_parameter("w52", w52);
    photon->set_shower_shape_parameter("w72", w72);
    photon->set_shower_shape_parameter("cluster_eta", cluster_eta);
    photon->set_shower_shape_parameter("cluster_phi", cluster_phi);

    // HCAL info
    std::vector<int> ihcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomIH, m_ihcal_tower_container, 0.0, true);
    std::vector<int> ohcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomOH, m_ohcal_tower_container, 0.0, false);

    float ihcal_et = 0;
    float ohcal_et = 0;
    float ihcal_et22 = 0;
    float ohcal_et22 = 0;
    float ihcal_et33 = 0;
    float ohcal_et33 = 0;

    int ihcal_ieta = ihcal_tower[0];
    int ihcal_iphi = ihcal_tower[1];
    float ihcalEt33[3][3] = {{0.0f}};

    int ohcal_ieta = ohcal_tower[0];
    int ohcal_iphi = ohcal_tower[1];
    float ohcalEt33[3][3] = {{0.0f}};

    for (int ieta_h = ihcal_ieta - 1; ieta_h <= ihcal_ieta + 1; ieta_h++)
    {
        for (int iphi_h = ihcal_iphi - 1; iphi_h <= ihcal_iphi + 1; iphi_h++)
        {
            int temp_ieta = ieta_h;
            int temp_iphi = iphi_h;
            shift_tower_index(temp_ieta, temp_iphi, 24, 64);
            if (temp_ieta < 0) continue;
            
            unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
            TowerInfo *towerinfo = m_ihcal_tower_container->get_tower_at_key(towerinfokey);
            if (towerinfo && towerinfo->get_isGood())
            {
                const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, temp_ieta, temp_iphi);
                RawTowerGeom *tower_geom = m_geomIH->get_tower_geometry(key);
                if(tower_geom)
                {
                    float energy = towerinfo->get_energy();
                    float eta = getTowerEta(tower_geom, 0, 0, 0);
                    float sintheta = 1.0 / cosh(eta);
                    float Et = energy * sintheta;
                    ihcalEt33[ieta_h - ihcal_ieta + 1][iphi_h - ihcal_iphi + 1] = Et;
                }
            }
        }
    }

    for (int ieta_h = ohcal_ieta - 1; ieta_h <= ohcal_ieta + 1; ieta_h++)
    {
        for (int iphi_h = ohcal_iphi - 1; iphi_h <= ohcal_iphi + 1; iphi_h++)
        {
            int temp_ieta = ieta_h;
            int temp_iphi = iphi_h;
            shift_tower_index(temp_ieta, temp_iphi, 24, 64);
            if (temp_ieta < 0) continue;

            unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
            TowerInfo *towerinfo = m_ohcal_tower_container->get_tower_at_key(towerinfokey);
            if (towerinfo && towerinfo->get_isGood())
            {
                 const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, temp_ieta, temp_iphi);
                 RawTowerGeom *tower_geom = m_geomOH->get_tower_geometry(key);
                 if(tower_geom)
                 {
                    float energy = towerinfo->get_energy();
                    float eta = getTowerEta(tower_geom, 0, 0, 0);
                    float sintheta = 1.0 / cosh(eta);
                    float Et = energy * sintheta;
                    ohcalEt33[ieta_h - ohcal_ieta + 1][iphi_h - ohcal_iphi + 1] = Et;
                 }
            }
        }
    }
    
    ihcal_et = ihcalEt33[1][1];
    ohcal_et = ohcalEt33[1][1];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            ihcal_et33 += ihcalEt33[i][j];
            ohcal_et33 += ohcalEt33[i][j];
            if (i == 1 || j == 1 + ihcal_tower[2])
            {
                if (j == 1 || i == 1 + ihcal_tower[3])
                {
                    ihcal_et22 += ihcalEt33[i][j];
                }
            }
            if (i == 1 || j == 1 + ohcal_tower[2])
            {
                if (j == 1 || i == 1 + ohcal_tower[3])
                {
                    ohcal_et22 += ohcalEt33[i][j];
                }
            }
        }
    }
    
    photon->set_shower_shape_parameter("ihcal_et", ihcal_et);
    photon->set_shower_shape_parameter("ohcal_et", ohcal_et);
    photon->set_shower_shape_parameter("ihcal_et22", ihcal_et22);
    photon->set_shower_shape_parameter("ohcal_et22", ohcal_et22);
    photon->set_shower_shape_parameter("ihcal_et33", ihcal_et33);
    photon->set_shower_shape_parameter("ohcal_et33", ohcal_et33);
    photon->set_shower_shape_parameter("ihcal_ieta", ihcal_ieta);
    photon->set_shower_shape_parameter("ihcal_iphi", ihcal_iphi);
    photon->set_shower_shape_parameter("ohcal_ieta", ohcal_ieta);
    photon->set_shower_shape_parameter("ohcal_iphi", ohcal_iphi);
}

double PhotonClusterBuilder::getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz)
{
  if (!tower_geom) return -9999;
  if (vx == 0 && vy == 0 && vz == 0)
  {
    return tower_geom->get_eta();
  }
  else
  {
    double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
    double theta = atan2(radius, tower_geom->get_center_z() - vz);
    return -log(tan(theta / 2.));
  }
}

std::vector<int> PhotonClusterBuilder::find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer *geom, TowerInfoContainer *towerContainer, float vertex_z, bool isihcal)
{
  int matchedieta = -1;
  int matchediphi = -1;
  double matchedeta = -999;
  double matchedphi = -999;
  
  if (!geom || !towerContainer) return {-1, -1, 0, 0};

  unsigned int ntowers = towerContainer->size();
  float minR = 999;

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo *tower = towerContainer->get_tower_at_channel(channel);
    if (!tower) continue;
    
    unsigned int towerkey = towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);
    
    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(isihcal ? RawTowerDefs::CalorimeterId::HCALIN : RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
    RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom) continue;

    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double dR_val = deltaR(eta, this_eta, phi, this_phi);
    if (dR_val < minR)
    {
      minR = dR_val;
      matchedieta = ieta;
      matchediphi = iphi;
      matchedeta = this_eta;
      matchedphi = this_phi;
    }
  }
  
  float deta = eta - matchedeta;
  float dphi_val = phi - matchedphi;
  if (dphi_val > M_PI) dphi_val -= 2 * M_PI;
  if (dphi_val < -M_PI) dphi_val += 2 * M_PI;
  
  int dphisign = (dphi_val > 0) ? 1 : -1;
  int detasign = (deta > 0) ? 1 : -1;

  return {matchedieta, matchediphi, detasign, dphisign};
}

double PhotonClusterBuilder::deltaR(double eta1, double phi1, double eta2, double phi2)
{
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi <= -M_PI) dphi += 2 * M_PI;
    return sqrt(pow(eta1 - eta2, 2) + pow(dphi, 2));
}
