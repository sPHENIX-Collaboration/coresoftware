#include "CaloVtxAlgoCaloZ.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <array>

CaloVtxAlgoCaloZ::CaloVtxAlgoCaloZ()
{

}

int CaloVtxAlgoCaloZ::Init(PHCompositeNode * /*topNode*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloVtxAlgoCaloZ::CalculateVertex(PHCompositeNode *topNode, float &zvtx)
{
  zvtx = std::numeric_limits<float>::quiet_NaN();


  TowerInfoContainer *emcal_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if (!emcal_towers || !hcalin_towers || !hcalout_towers)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

  if (!tower_geomEM || !tower_geomIH || !tower_geomOH)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  int size;

  float average_z[3]{0};
  float total_E[3]{0};

  if (emcal_towers)
  {
    size = emcal_towers->size();  // online towers should be the same!
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo *_tower = emcal_towers->get_tower_at_channel(channel);
      short good = (_tower->get_isGood() ? 1 : 0);
      if (!good)
      {
        continue;
      }

      float energy = _tower->get_energy();
      if (energy < m_energy_cut)
      {
        continue;
      }

      // float time = _tower->get_time_float();

      unsigned int towerkey = emcal_towers->encode_key(channel);
      int ieta = emcal_towers->getTowerEtaBin(towerkey);
      int iphi = emcal_towers->getTowerPhiBin(towerkey);

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
      float tower_z = tower_geomEM->get_tower_geometry(key)->get_center_z();
      average_z[0] += tower_z * energy;
      total_E[0] += energy;
    }
  }

  if (hcalin_towers)
  {
    size = hcalin_towers->size();  // online towers should be the same!
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo *_tower = hcalin_towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      if (energy < m_energy_cut)
      {
        continue;
      }
      // float time = _tower->get_time_float();
      short good = (_tower->get_isGood() ? 1 : 0);
      if (!good)
      {
        continue;
      }

      unsigned int towerkey = hcalin_towers->encode_key(channel);
      int ieta = hcalin_towers->getTowerEtaBin(towerkey);
      int iphi = hcalin_towers->getTowerPhiBin(towerkey);
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
      float tower_z = tower_geomIH->get_tower_geometry(key)->get_center_z();
      average_z[1] += tower_z * energy;
      total_E[1] += energy;
    }
  }
  if (hcalout_towers)
  {
    size = hcalout_towers->size();  // online towers should be the same!
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo *_tower = hcalout_towers->get_tower_at_channel(channel);
      float energy = _tower->get_energy();
      if (energy < m_energy_cut)
      {
        continue;
      }
      // float time = _tower->get_time_float();
      unsigned int towerkey = hcalout_towers->encode_key(channel);
      int ieta = hcalout_towers->getTowerEtaBin(towerkey);
      int iphi = hcalout_towers->getTowerPhiBin(towerkey);
      short good = (_tower->get_isGood() ? 1 : 0);

      if (!good)
      {
        continue;
      }

      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);

      float tower_z = tower_geomOH->get_tower_geometry(key)->get_center_z();
      average_z[2] += tower_z * energy;
      total_E[2] += energy;
    }
  }

  double b_calo_vertex_z = (average_z[0] + average_z[1] + average_z[2]) / (total_E[0] + total_E[1] + total_E[2]);

  zvtx = b_calo_vertex_z;

  return Fun4AllReturnCodes::EVENT_OK;
}
