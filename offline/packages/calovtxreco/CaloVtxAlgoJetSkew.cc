#include "CaloVtxAlgoJetSkew.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <array>

CaloVtxAlgoJetSkew::CaloVtxAlgoJetSkew()
{
}

int CaloVtxAlgoJetSkew::Init(PHCompositeNode *topNode)
{
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if(!tower_geomEM || !tower_geomOH)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  m_radius_EM = tower_geomEM->get_radius();
  m_radius_OH = tower_geomOH->get_radius();

  if(std::isnan(m_radius_EM) || std::isnan(m_radius_OH))
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloVtxAlgoJetSkew::CalculateVertex(PHCompositeNode *topNode, float &zvtx)
{

  zvtx = std::numeric_limits<float>::quiet_NaN();

  JetContainer *jetcon = findNode::getClass<JetContainer>(topNode, m_jetnodename);

  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  const int nz = 601;
  const int njet = 2;
  Jet *jets[njet];
  float jpt[njet] = {0};
  float jemsum[njet] = {0};
  float johsum[njet] = {0};
  float jemeta[njet] = {0};
  float joheta[njet] = {0};

  if (jetcon)
  {
    int tocheck = jetcon->size();
    for (int i = 0; i < tocheck; ++i)
    {
      Jet *jet = jetcon->get_jet(i);
      if (jet)
      {
        float pt = jet->get_pt();
        if (pt < m_jet_threshold)
        {
          continue;
        }
        if (pt > jpt[0])
        {
          jpt[1] = jpt[0];
          jets[1] = jets[0];
          jets[0] = jet;
          jpt[0] = pt;
        }
        else if (pt > jpt[1])
        {
          jets[1] = jet;
          jpt[1] = pt;
        }
      }
    }
  }
  else
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float metric = std::numeric_limits<float>::max();
  for (int i = 0; i < nz; ++i)
  {
    float testz = -300 + i;
    float testmetric = 0;
    for (int j = 0; j < njet; ++j)
    {
      if (jpt[j] == 0)
      {
        continue;
      }
      jemsum[j] = 0;
      johsum[j] = 0;
      jemeta[j] = 0;
      joheta[j] = 0;

      for (auto comp : jets[j]->get_comp_vec())
	{
	  if (comp.first == 5 || comp.first == 26)
	    {
	      continue;
	    }
	  unsigned int channel = comp.second;
	  if (comp.first == 7 || comp.first == 27)
	    {
	      TowerInfo *tower = towers[2]->get_tower_at_channel(channel);
	      if (tower->get_energy() < 0.1)
		{
		  continue;
		}
	      johsum[j] += tower->get_energy();
	      float neweta = new_eta(channel, towers[2], geom[2], RawTowerDefs::CalorimeterId::HCALOUT, testz);
	      joheta[j] += neweta * tower->get_energy();
	    }
	  if (comp.first == 13 || comp.first == 28 || comp.first == 25)
	    {
	      TowerInfo *tower = towers[0]->get_tower_at_channel(channel);
	      if (tower->get_energy() < 0.1)
		{
		  continue;
		}
	      jemsum[j] += tower->get_energy();
	      float neweta = new_eta(channel, towers[0], geom[1], RawTowerDefs::CalorimeterId::HCALIN, testz);
	      jemeta[j] += neweta * tower->get_energy();
	    }
	}
    
      jemeta[j] /= jemsum[j];
      joheta[j] /= johsum[j];
      if (!std::isnan(jemeta[j]) && !std::isnan(joheta[j]))
	{
	  testmetric += pow(jemeta[j] - joheta[j], 2);
	}
    }
    if (testmetric < metric && testmetric != 0)
      {
	metric = testmetric;
	zvtx = testz;
      }
  }
  if (fabs(zvtx) >= 305)
    {
      zvtx = std::numeric_limits<float>::quiet_NaN();
    }
  else
    {
      zvtx *= m_calib_factor;  // calibration factor from simulation
    }
  
  return Fun4AllReturnCodes::EVENT_OK;

}

float CaloVtxAlgoJetSkew::new_eta(int channel, TowerInfoContainer *towerset, RawTowerGeomContainer *geom, RawTowerDefs::CalorimeterId caloID, float testz)
{
  int key = towerset->encode_key(channel);
  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towerset->getTowerEtaBin(key), towerset->getTowerPhiBin(key));

  RawTowerGeom *tower_geom = geom->get_tower_geometry(geomkey);
  float oldeta = tower_geom->get_eta();

  float radius = (caloID == RawTowerDefs::CalorimeterId::HCALIN ? m_radius_EM : m_radius_OH);
  float towerz = radius / (tanf(2 * std::atanf(std::exp(oldeta))));
  float newz = towerz + testz;
  float newTheta = std::atan2(radius, newz);
  float neweta = -log(tan(0.5 * newTheta));

  return neweta;
}
