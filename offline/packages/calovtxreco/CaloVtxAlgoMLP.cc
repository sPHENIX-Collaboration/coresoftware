#include "CaloVtxAlgoMLP.h"
#include "VertexMLP.h"

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

CaloVtxAlgoMLP::CaloVtxAlgoMLP()
{
}

int CaloVtxAlgoMLP::Init(PHCompositeNode *topNode)
{
  if (!VertexMLP::Load(m_weightsFile))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
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

int CaloVtxAlgoMLP::CalculateVertex(PHCompositeNode *topNode, float &zvtx)
{
  TowerInfoContainer *emcalre_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  TowerInfoContainer *hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if (!emcalre_towers || !hcalin_towers || !hcalout_towers)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

  if (!tower_geomEM || !tower_geomIH || !tower_geomOH)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // TODO: fill `features` from topNode in the exact order documented in
  // VertexMLP.h (emcal/ohcal lead+sublead flags, zmean/zsig/zskew,
  // energy, exj). This is the same per-event extraction CaloVtxReco used
  // to do inline before the NN call -- move that block here unchanged.

  JetContainer *jetscon = findNode::getClass<JetContainer>(topNode, m_jet_node.c_str());

  float dijet_pt[2]={0};
  float dijet_eta[2]={0};
  float dijet_phi[2]={0};
  float dijet_E[2]={0};
  float dijet_em_frac[2]={0};
  float dijet_em_mean[2]={0};
  float dijet_em_skew[2]={0};
  float dijet_em_sig[2]={0};
  float dijet_oh_frac[2]={0};
  float dijet_oh_mean[2]={0};
  float dijet_oh_skew[2]={0};
  float dijet_oh_sig[2]={0};

  if (!jetscon)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  for (auto jet : *jetscon)
    {
      double jet_E = 0;
      double jet_emcal = 0;
      double jet_ohcal = 0;
      double jetpt = jet->get_pt();
      double jetphi = jet->get_phi();
      double jeteta= jet->get_eta();

      if (jetpt < 0.5) continue;

      double em_S0 = 0.0;
      double em_S1 = 0.0;
      double em_S2 = 0.0;
      double em_S3 = 0.0;
      double oh_S0 = 0.0;
      double oh_S1 = 0.0;
      double oh_S2 = 0.0;
      double oh_S3 = 0.0;


      int itower = 0;
      for (auto comp : jet->get_comp_vec())
	{
	  unsigned int channel = comp.second;
	  TowerInfo *tower;
	  float tower_e = 0;
	  if (comp.first == 26 || comp.first == 30)
	    {
	      tower = hcalin_towers->get_tower_at_channel(channel);
		      
	      if (!tower || !tower_geomIH)
		{
		  continue;
		}
		      
	      tower_e = tower->get_energy();
	      if (tower_e < 0.005) continue;
	      jet_E += tower_e;
	    }
	  else if (comp.first == 27 || comp.first == 31)
	    {
	      // OHCAL
	      tower = hcalout_towers->get_tower_at_channel(channel);
	      if (!tower || !tower_geomOH)
		{
		  continue;
		}
	      unsigned int towerkey = hcalout_towers->encode_key(channel);
	      int ieta = hcalout_towers->getTowerEtaBin(towerkey);
	      int iphi = hcalout_towers->getTowerPhiBin(towerkey);
	      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
	      float tower_eta = tower_geomOH->get_tower_geometry(key)->get_eta();
	      float tower_r = m_radius_OH;
	      float tower_z = tower_r*sinh(tower_eta);
		      
	      tower_e = tower->get_energy();
	      if (tower_e < 0.035) continue;
	      jet_ohcal += tower_e;
	      oh_S0 += tower_e;
	      oh_S1 += tower_e*tower_z;
	      oh_S2 += tower_e*tower_z*tower_z;
	      oh_S3 += tower_e*tower_z*tower_z*tower_z;
	      jet_E += tower_e;
	    }	       
	  else if (comp.first == 28 || comp.first == 29)
	    {  // EMCAL
	      tower = emcalre_towers->get_tower_at_channel(channel);

	      if (!tower || !tower_geomOH)
		{
		  continue;
		}
	      unsigned int towerkey = emcalre_towers->encode_key(channel);
	      int ieta = emcalre_towers->getTowerEtaBin(towerkey);
	      int iphi = emcalre_towers->getTowerPhiBin(towerkey);
	      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
	      float tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();
	      float tower_r = m_radius_EM;
	      float tower_z = tower_r*sinh(tower_eta);
		      
	      tower_e = tower->get_energy();
	      if (tower_e < 0.068) continue;
	      jet_emcal += tower_e;
	      jet_E += tower_e;
	      em_S0 += tower_e;
	      em_S1 += tower_e*tower_z;
	      em_S2 += tower_e*tower_z*tower_z;
	      em_S3 += tower_e*tower_z*tower_z*tower_z;
	    }

	  itower++;	      
	}

      if (itower == 0) continue;

      double em_mean = em_S1/em_S0;

      double em_variance = em_S2/em_S0 - em_mean*em_mean;
      double em_sigma = std::sqrt(std::max(0.0, em_variance));

      double em_mu3 =
	em_S3/em_S0
	- 3.0*em_mean*(em_S2/em_S0)
	+ 2.0*em_mean*em_mean*em_mean;

      double em_skew = (em_sigma > 0.0)
	? em_mu3/std::pow(em_sigma,3)
	: 0.0;

      double oh_mean = oh_S1/oh_S0;

      double oh_variance = oh_S2/oh_S0 - oh_mean*oh_mean;
      double oh_sigma = std::sqrt(std::max(0.0, oh_variance));

      double oh_mu3 =
	oh_S3/oh_S0
	- 3.0*oh_mean*(oh_S2/oh_S0)
	+ 2.0*oh_mean*oh_mean*oh_mean;

      double oh_skew = (oh_sigma > 0.0)
	? oh_mu3/std::pow(oh_sigma,3)
	: 0.0;
	      
      jet_emcal /= jet_E;
      jet_ohcal /= jet_E;

      if (jetpt > dijet_pt[0])
	{
	  dijet_pt[1]=dijet_pt[0];
	  dijet_eta[1]= dijet_eta[0];
	  dijet_phi[1]= dijet_phi[0];
	  dijet_E[1]= dijet_E[0];
	  dijet_em_frac[1]= dijet_em_frac[0];
	  dijet_em_mean[1]= dijet_em_mean[0];
	  dijet_em_skew[1]= dijet_em_skew[0];
	  dijet_em_sig[1]= dijet_em_sig[0];
	  dijet_oh_frac[1]= dijet_oh_frac[0];
	  dijet_oh_mean[1]= dijet_oh_mean[0];
	  dijet_oh_skew[1]= dijet_oh_skew[0];
	  dijet_oh_sig[1]= dijet_oh_sig[0];

	  dijet_pt[0]=jetpt;
	  dijet_eta[0]= jeteta;
	  dijet_phi[0]= jetphi;
	  dijet_E[0]= jet_E;
	  dijet_em_frac[0]= jet_emcal;
	  dijet_em_mean[0]= em_mean;
	  dijet_em_skew[0]= em_skew;
	  dijet_em_sig[0]= em_sigma;
	  dijet_oh_frac[0]= jet_ohcal;
	  dijet_oh_mean[0]= oh_mean;
	  dijet_oh_skew[0]= oh_skew;
	  dijet_oh_sig[0]= oh_sigma;

	}
      else if (jetpt > dijet_pt[1])
	{
	  dijet_pt[1]=jetpt;
	  dijet_eta[1]= jeteta;
	  dijet_phi[1]= jetphi;
	  dijet_E[1]= jet_E;
	  dijet_em_frac[1]= jet_emcal;
	  dijet_em_mean[1]= em_mean;
	  dijet_em_skew[1]= em_skew;
	  dijet_em_sig[1]= em_sigma;
	  dijet_oh_frac[1]= jet_ohcal;
	  dijet_oh_mean[1]= oh_mean;
	  dijet_oh_skew[1]= oh_skew;
	  dijet_oh_sig[1]= oh_sigma;
	}
    }

  double emcal_on[2] = {1};
  double ohcal_on[2] = {1};

	  
  for (int i = 0 ; i < 2; i++)
    {
      if (std::isnan(dijet_em_frac[i]) ||
	  std::isnan(dijet_em_mean[i]) ||
	  std::isnan(dijet_em_skew[i]) ||
	  std::isnan(dijet_em_sig[i]))
	{
	  emcal_on[i] = 0;
	}
      if (std::isnan(dijet_oh_frac[i]) ||
	  std::isnan(dijet_oh_mean[i]) ||
	  std::isnan(dijet_oh_skew[i]) ||
	  std::isnan(dijet_oh_sig[i]))
	{
	  ohcal_on[i] = 0;
	}
    }
		     
  float emcal_z_moments[2][4] = {0};
  float ohcal_z_moments[2][4] = {0};
  for (int i = 0; i < 2; i++)
    {
      if (emcal_on[i])
	{
	  emcal_z_moments[i][0] = dijet_em_frac[i];
	  emcal_z_moments[i][1] = dijet_em_mean[i];
	  emcal_z_moments[i][2] = dijet_em_skew[i];
	  emcal_z_moments[i][3] = dijet_em_sig[i];
	}
      if (ohcal_on[i])
	{
	  ohcal_z_moments[i][0] = dijet_oh_frac[i];
	  ohcal_z_moments[i][1] = dijet_oh_mean[i];
	  ohcal_z_moments[i][2] = dijet_oh_skew[i];
	  ohcal_z_moments[i][3] = dijet_oh_sig[i];
	}
    }

  float exj = dijet_E[1]/dijet_E[0];
  
  std::array<double, VertexMLP::kNFeatures> features = {
    emcal_on[0],
    ohcal_on[0],
    emcal_on[1],
    ohcal_on[1],
    emcal_z_moments[0][0],
    emcal_z_moments[0][1],
    emcal_z_moments[0][2],
    emcal_z_moments[0][3],
    ohcal_z_moments[0][0],
    ohcal_z_moments[0][1],
    ohcal_z_moments[0][2],
    ohcal_z_moments[0][3],
    emcal_z_moments[1][0],
    emcal_z_moments[1][1],
    emcal_z_moments[1][2],
    emcal_z_moments[1][3],
    ohcal_z_moments[1][0],
    ohcal_z_moments[1][1],
    ohcal_z_moments[1][2],
    ohcal_z_moments[1][3],
    exj};

  
  zvtx = static_cast<float>(VertexMLP::PredictVertexZ(features));
  return Fun4AllReturnCodes::EVENT_OK;
}
