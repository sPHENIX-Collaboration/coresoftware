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

int CaloVtxAlgoMLP::Init(PHCompositeNode *topNode)
{
  if (!VertexMLP::Load(m_weightsFile))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  RawTowerGeomContainer *tower_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *tower_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  RawTowerGeomContainer *tower_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if(!tower_geomEM || !tower_geomOH)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
  unsigned int emkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 48, 128);
  m_radius_EM = tower_geomEM->get_tower_geometry(emkey)->get_center_radius();
  unsigned int ihkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, 12, 0);
  m_radius_IH = tower_geomIH->get_tower_geometry(ihkey)->get_center_radius();
  unsigned int ohkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, 12, 0);
  m_radius_OH = tower_geomOH->get_tower_geometry(ohkey)->get_center_radius();
  
  // m_radius_EM = tower_geomEM->get_radius();
  // m_radius_IH = tower_geomIH->get_radius();
  // m_radius_OH = tower_geomOH->get_radius();

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


  std::vector<float> *b_reco_jet_pt = new std::vector<float>();
  std::vector<float> *b_reco_jet_eta = new std::vector<float>();
  std::vector<float> *b_reco_jet_phi = new std::vector<float>();
  std::vector<float> *b_reco_jet_emcal_eta = new std::vector<float>();
  std::vector<float> *b_reco_jet_ihcal_eta = new std::vector<float>();
  std::vector<float> *b_reco_jet_ohcal_eta = new std::vector<float>();
  std::vector<float> *b_reco_jet_emcal_zmean = new std::vector<float>();
  std::vector<float> *b_reco_jet_ohcal_zmean = new std::vector<float>();
  std::vector<float> *b_reco_jet_emcal_zsig = new std::vector<float>();
  std::vector<float> *b_reco_jet_ohcal_zsig = new std::vector<float>();
  std::vector<float> *b_reco_jet_emcal_zskew = new std::vector<float>();
  std::vector<float> *b_reco_jet_ohcal_zskew = new std::vector<float>();
  std::vector<float> *b_reco_jet_emcal = new std::vector<float>();
  std::vector<float> *b_reco_jet_ohcal = new std::vector<float>();
  std::vector<float> *b_reco_jet_ihcal = new std::vector<float>();

  JetContainer *jetscon = findNode::getClass<JetContainer>(topNode, m_jet_node.c_str());

  int ijet = 0;
  double dijet_pt[2] = {0};
  int dijet_index[2] = {0};

  if (jetscon)
    {
      for (auto jet : *jetscon)
	{
	  double jet_E = 0;
	  double jet_emcal = 0;
	  double jet_ihcal = 0;
	  double jet_ohcal = 0;
	  double jet_emcal_eta = 0;
	  double jet_ihcal_eta = 0;
	  double jet_ohcal_eta = 0;
	  double jetpt = jet->get_pt();
	  if (jetpt < 2) continue;

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
		{  // IHcal
		  tower = hcalin_towers->get_tower_at_channel(channel);
		      
		  if (!tower || !tower_geomIH)
		    {
		      continue;
		    }
		      
		  unsigned int towerkey = hcalin_towers->encode_key(channel);
		  int ieta = hcalin_towers->getTowerEtaBin(towerkey);
		  int iphi = hcalin_towers->getTowerPhiBin(towerkey);
		  const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
		  float tower_eta = tower_geomIH->get_tower_geometry(key)->get_eta();
		  float tower_r = m_radius_IH;
		  float tower_z = tower_r*sinh(tower_eta);
		  float new_tower_z = tower_z;// - jet_vertex;
		  float new_tower_eta = asinh(new_tower_z/tower_r);
		      
		  tower_e = tower->get_energy();
		  if (tower_e < 0.005) continue;
		  jet_ihcal_eta += new_tower_eta*tower_e;
		  jet_ihcal += tower_e;
		  jet_E += tower_e;
		}
	      else if (comp.first == 27 || comp.first == 31)
		{  // OHCAL
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
		  float new_tower_z = tower_z;// - jet_vertex;
		  float new_tower_eta = asinh(new_tower_z/tower_r);
		      
		  tower_e = tower->get_energy();
		  if (tower_e < 0.035) continue;
		  jet_ohcal_eta += new_tower_eta*tower_e;
		  jet_ohcal += tower_e;
		  oh_S0 += tower_e;
		  oh_S1 += tower_e*tower_z;
		  oh_S2 += tower_e*tower_z*tower_z;
		  oh_S3 += tower_e*tower_z*tower_z*tower_z;

		  //tower_e = tower->get_energy();
		  //jet_ohcal += tower_e;
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
		  float new_tower_z = tower_z;// - jet_vertex;
		  float new_tower_eta = asinh(new_tower_z/tower_r);
		      
		  tower_e = tower->get_energy();
		  if (tower_e < 0.068) continue;
		  jet_emcal_eta += new_tower_eta*tower_e;
		  jet_emcal += tower_e;
		  jet_E += tower_e;
		  em_S0 += tower_e;
		  em_S1 += tower_e*tower_z;
		  em_S2 += tower_e*tower_z*tower_z;
		  em_S3 += tower_e*tower_z*tower_z*tower_z;

		  // if (print)
		  // 	{
		  // 	  std::cout << itower << " | e:" << tower_e <<  " / eta: " << tower_eta <<  " / " << new_tower_eta << std::endl;
		  // 	}
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
	      
	  jet_emcal_eta /= jet_emcal;
	  jet_ihcal_eta /= jet_ihcal;
	  jet_ohcal_eta /= jet_ohcal;
	  jet_emcal /= jet_E;
	  jet_ihcal /= jet_E;
	  jet_ohcal /= jet_E;
	  

	  b_reco_jet_pt->push_back(jet->get_pt());
	  b_reco_jet_eta->push_back(jet->get_eta());
	  b_reco_jet_phi->push_back(jet->get_phi());

	  b_reco_jet_emcal->push_back(jet_emcal);
	  b_reco_jet_ohcal->push_back(jet_ohcal);
	  b_reco_jet_ihcal->push_back(jet_ihcal);

	  b_reco_jet_emcal_zmean->push_back(em_mean);
	  b_reco_jet_ohcal_zmean->push_back(oh_mean);
	  b_reco_jet_emcal_zsig->push_back(em_sigma);
	  b_reco_jet_ohcal_zsig->push_back(oh_sigma);
	  b_reco_jet_emcal_zskew->push_back(em_skew);
	  b_reco_jet_ohcal_zskew->push_back(oh_skew);


	  b_reco_jet_emcal_eta->push_back(jet_emcal_eta);
	  b_reco_jet_ohcal_eta->push_back(jet_ohcal_eta);
	  b_reco_jet_ihcal_eta->push_back(jet_ihcal_eta);


	  if (jetpt > dijet_pt[0])
	    {
	      dijet_pt[1] = dijet_pt[0];
	      dijet_index[1] = dijet_index[0];
	      dijet_pt[0] = jetpt;
	      dijet_index[0] = ijet;
	    }
	  else if (jetpt > dijet_pt[1])
	    {
	      dijet_pt[1] = jetpt;
	      dijet_index[1] = ijet;
	    }
	  ijet++;
	}

    }

  // TODO: fill `features` from topNode in the exact order documented in
  // VertexMLP.h (emcal/ohcal lead+sublead flags, zmean/zsig/zskew,
  // energy, exj). This is the same per-event extraction CaloVtxReco used
  // to do inline before the NN call -- move that block here unchanged.
  float e1 = dijet_pt[0] * cosh(b_reco_jet_eta->at(dijet_index[0]));
  float e2 = dijet_pt[1] * cosh(b_reco_jet_eta->at(dijet_index[1]));

  double exj = e2/e1;
  double emcal_lead_on = 1;
  double emcal_sublead_on = 1;
  double ohcal_lead_on = 1;
  double ohcal_sublead_on = 1;
	  
  if (std::isnan(b_reco_jet_emcal_zmean->at(dijet_index[0]))) emcal_lead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal_zsig->at(dijet_index[0]))) emcal_lead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal_zskew->at(dijet_index[0]))) emcal_lead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal->at(dijet_index[0]))) emcal_lead_on = 0;//<","
  if (std::isnan(b_reco_jet_ohcal_zmean->at(dijet_index[0]))) ohcal_lead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal_zsig->at(dijet_index[0]))) ohcal_lead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal_zskew->at(dijet_index[0]))) ohcal_lead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal->at(dijet_index[0]))) ohcal_lead_on = 0;//<","
  if (std::isnan(b_reco_jet_emcal_zmean->at(dijet_index[1]))) emcal_sublead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal_zsig->at(dijet_index[1]))) emcal_sublead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal_zskew->at(dijet_index[1]))) emcal_sublead_on = 0;//maxz_EM <<","
  if (std::isnan(b_reco_jet_emcal->at(dijet_index[1]))) emcal_sublead_on = 0;//<","
  if (std::isnan(b_reco_jet_ohcal_zmean->at(dijet_index[1]))) ohcal_sublead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal_zsig->at(dijet_index[1]))) ohcal_sublead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal_zskew->at(dijet_index[1]))) ohcal_sublead_on = 0;//maxz_OH <<","
  if (std::isnan(b_reco_jet_ohcal->at(dijet_index[1]))) ohcal_sublead_on = 0;//<","

  float emcal_lead_z_moments[4] = {b_reco_jet_emcal_zmean->at(dijet_index[0]),
    b_reco_jet_emcal_zsig->at(dijet_index[0]),
    b_reco_jet_emcal_zskew->at(dijet_index[0]),
    b_reco_jet_emcal->at(dijet_index[0])};
  float emcal_sublead_z_moments[4] = {b_reco_jet_emcal_zmean->at(dijet_index[1]),
    b_reco_jet_emcal_zsig->at(dijet_index[1]),
    b_reco_jet_emcal_zskew->at(dijet_index[1]),
    b_reco_jet_emcal->at(dijet_index[1])};
  float ohcal_lead_z_moments[4] = {b_reco_jet_ohcal_zmean->at(dijet_index[0]),
    b_reco_jet_ohcal_zsig->at(dijet_index[0]),
    b_reco_jet_ohcal_zskew->at(dijet_index[0]),
    b_reco_jet_ohcal->at(dijet_index[0])};
  float ohcal_sublead_z_moments[4] = {b_reco_jet_ohcal_zmean->at(dijet_index[1]),
    b_reco_jet_ohcal_zsig->at(dijet_index[1]),
    b_reco_jet_ohcal_zskew->at(dijet_index[1]),
    b_reco_jet_ohcal->at(dijet_index[1])};

  if (!emcal_lead_on)
    {
      for (int i = 0; i < 4; i++)
	{
	  emcal_lead_z_moments[i] = 0;
	}
    }
  if (!emcal_sublead_on)
    {
      for (int i = 0; i < 4; i++)
	{
	  emcal_sublead_z_moments[i] = 0;
	}
    }
  if (!ohcal_lead_on)
    {
      for (int i = 0; i < 4; i++)
	{
	  ohcal_lead_z_moments[i] = 0;
	}
    }
  if (!ohcal_sublead_on)
    {
      for (int i = 0; i < 4; i++)
	{
	  ohcal_sublead_z_moments[i] = 0;
	}
    }
  // {
  //   std::cout << emcal_lead_on << ","
  // 	    << ohcal_lead_on << ","
  // 	    << emcal_sublead_on << ","
  // 	    << ohcal_sublead_on << ","
  // 	    << emcal_lead_z_moments[0] << ","
  // 	    << emcal_lead_z_moments[1] << ","
  // 	    << emcal_lead_z_moments[2] << ","
  // 	    << emcal_lead_z_moments[3] << ","
  // 	    << ohcal_lead_z_moments[0] << ","
  // 	    << ohcal_lead_z_moments[1] << ","
  // 	    << ohcal_lead_z_moments[2] << ","
  // 	    << ohcal_lead_z_moments[3] << ","
  // 	    << emcal_sublead_z_moments[0] << ","
  // 	    << emcal_sublead_z_moments[1] << ","
  // 	    << emcal_sublead_z_moments[2] << ","
  // 	    << emcal_sublead_z_moments[3] << ","
  // 	    << ohcal_sublead_z_moments[0] << ","
  // 	    << ohcal_sublead_z_moments[1] << ","
  // 	    << ohcal_sublead_z_moments[2] << ","
  // 	    << ohcal_sublead_z_moments[3] << ","
  // 	    << exj << std::endl;
  // }

  std::array<double, VertexMLP::kNFeatures> features = {emcal_lead_on,
    ohcal_lead_on,
    emcal_sublead_on,
    ohcal_sublead_on,
    emcal_lead_z_moments[0],
    emcal_lead_z_moments[1],
    emcal_lead_z_moments[2],
    emcal_lead_z_moments[3],
    ohcal_lead_z_moments[0],
    ohcal_lead_z_moments[1],
    ohcal_lead_z_moments[2],
    ohcal_lead_z_moments[3],
    emcal_sublead_z_moments[0],
    emcal_sublead_z_moments[1],
    emcal_sublead_z_moments[2],
    emcal_sublead_z_moments[3],
    ohcal_sublead_z_moments[0],
    ohcal_sublead_z_moments[1],
    ohcal_sublead_z_moments[2],
    ohcal_sublead_z_moments[3],
    exj};

  
  zvtx = static_cast<float>(VertexMLP::PredictVertexZ(features));

  return Fun4AllReturnCodes::EVENT_OK;
}
