#include "ActsTrkFitAnalyzer.h"

/// General fun4all and subsysreco includes
#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TTree.h>


ActsTrkFitAnalyzer::ActsTrkFitAnalyzer(const std::string& name)
  : SubsysReco(name)
{
}

ActsTrkFitAnalyzer::~ActsTrkFitAnalyzer()
{
}

int ActsTrkFitAnalyzer::Init(PHCompositeNode *topNode)
{
  initializeTree();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int ActsTrkFitAnalyzer::process_event(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int ActsTrkFitAnalyzer::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

void ActsTrkFitAnalyzer::initializeTree()
{

  m_trackFile = new TFile(Name().c_str(),"RECREATE");
  
  m_trackTree = new TTree("tracktree","A tree with Acts KF track information");
  /*
  m_trackTree->Branch("event_nr", &m_eventNr);
  m_trackTree->Branch("traj_nr", &m_trajNr);
  m_trackTree->Branch("t_barcode", &m_t_barcode, "t_barcode/l");
  m_trackTree->Branch("t_charge", &m_t_charge);
  m_trackTree->Branch("t_time", &m_t_time);
  m_trackTree->Branch("t_vx", &m_t_vx);
  m_trackTree->Branch("t_vy", &m_t_vy);
  m_trackTree->Branch("t_vz", &m_t_vz);
  m_trackTree->Branch("t_px", &m_t_px);
  m_trackTree->Branch("t_py", &m_t_py);
  m_trackTree->Branch("t_pz", &m_t_pz);
  m_trackTree->Branch("t_theta", &m_t_theta);
  m_trackTree->Branch("t_phi", &m_t_phi);
  m_trackTree->Branch("t_eta", &m_t_eta);
  m_trackTree->Branch("t_pT", &m_t_pT);
  
  m_trackTree->Branch("t_x", &m_t_x);
  m_trackTree->Branch("t_y", &m_t_y);
  m_trackTree->Branch("t_z", &m_t_z);
  m_trackTree->Branch("t_r", &m_t_r);
  m_trackTree->Branch("t_dx", &m_t_dx);
  m_trackTree->Branch("t_dy", &m_t_dy);
  m_trackTree->Branch("t_dz", &m_t_dz);
  m_trackTree->Branch("t_eLOC0", &m_t_eLOC0);
  m_trackTree->Branch("t_eLOC1", &m_t_eLOC1);
  m_trackTree->Branch("t_ePHI", &m_t_ePHI);
  m_trackTree->Branch("t_eTHETA", &m_t_eTHETA);
  m_trackTree->Branch("t_eQOP", &m_t_eQOP);
  m_trackTree->Branch("t_eT", &m_t_eT);
  
  m_trackTree->Branch("nStates", &m_nStates);
  m_trackTree->Branch("nMeasurements", &m_nMeasurements);
  m_trackTree->Branch("volume_id", &m_volumeID);
  m_trackTree->Branch("layer_id", &m_layerID);
  m_trackTree->Branch("module_id", &m_moduleID);
  m_trackTree->Branch("l_x_hit", &m_lx_hit);
  m_trackTree->Branch("l_y_hit", &m_ly_hit);
  m_trackTree->Branch("g_x_hit", &m_x_hit);
  m_trackTree->Branch("g_y_hit", &m_y_hit);
  m_trackTree->Branch("g_z_hit", &m_z_hit);
  m_trackTree->Branch("res_x_hit", &m_res_x_hit);
  m_trackTree->Branch("res_y_hit", &m_res_y_hit);
  m_trackTree->Branch("err_x_hit", &m_err_x_hit);
  m_trackTree->Branch("err_y_hit", &m_err_y_hit);
  m_trackTree->Branch("pull_x_hit", &m_pull_x_hit);
  m_trackTree->Branch("pull_y_hit", &m_pull_y_hit);
  m_trackTree->Branch("dim_hit", &m_dim_hit);
  
  m_trackTree->Branch("hasFittedParams", &m_hasFittedParams);
  m_trackTree->Branch("eLOC0_fit", &m_eLOC0_fit);
  m_trackTree->Branch("eLOC1_fit", &m_eLOC1_fit);
  m_trackTree->Branch("ePHI_fit", &m_ePHI_fit);
  m_trackTree->Branch("eTHETA_fit", &m_eTHETA_fit);
  m_trackTree->Branch("eQOP_fit", &m_eQOP_fit);
  m_trackTree->Branch("eT_fit", &m_eT_fit);
  m_trackTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
  m_trackTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
  m_trackTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
  m_trackTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
  m_trackTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
  m_trackTree->Branch("err_eT_fit", &m_err_eT_fit);
  
  m_trackTree->Branch("nPredicted", &m_nPredicted);
  m_trackTree->Branch("predicted", &m_prt);
  m_trackTree->Branch("eLOC0_prt", &m_eLOC0_prt);
  m_trackTree->Branch("eLOC1_prt", &m_eLOC1_prt);
  m_trackTree->Branch("ePHI_prt", &m_ePHI_prt);
  m_trackTree->Branch("eTHETA_prt", &m_eTHETA_prt);
  m_trackTree->Branch("eQOP_prt", &m_eQOP_prt);
  m_trackTree->Branch("eT_prt", &m_eT_prt);
  m_trackTree->Branch("res_eLOC0_prt", &m_res_eLOC0_prt);
  m_trackTree->Branch("res_eLOC1_prt", &m_res_eLOC1_prt);
  m_trackTree->Branch("res_ePHI_prt", &m_res_ePHI_prt);
  m_trackTree->Branch("res_eTHETA_prt", &m_res_eTHETA_prt);
  m_trackTree->Branch("res_eQOP_prt", &m_res_eQOP_prt);
  m_trackTree->Branch("res_eT_prt", &m_res_eT_prt);
  m_trackTree->Branch("err_eLOC0_prt", &m_err_eLOC0_prt);
  m_trackTree->Branch("err_eLOC1_prt", &m_err_eLOC1_prt);
  m_trackTree->Branch("err_ePHI_prt", &m_err_ePHI_prt);
  m_trackTree->Branch("err_eTHETA_prt", &m_err_eTHETA_prt);
  m_trackTree->Branch("err_eQOP_prt", &m_err_eQOP_prt);
  m_trackTree->Branch("err_eT_prt", &m_err_eT_prt);
  m_trackTree->Branch("pull_eLOC0_prt", &m_pull_eLOC0_prt);
  m_trackTree->Branch("pull_eLOC1_prt", &m_pull_eLOC1_prt);
  m_trackTree->Branch("pull_ePHI_prt", &m_pull_ePHI_prt);
  m_trackTree->Branch("pull_eTHETA_prt", &m_pull_eTHETA_prt);
  m_trackTree->Branch("pull_eQOP_prt", &m_pull_eQOP_prt);
  m_trackTree->Branch("pull_eT_prt", &m_pull_eT_prt);
  m_trackTree->Branch("g_x_prt", &m_x_prt);
  m_trackTree->Branch("g_y_prt", &m_y_prt);
  m_trackTree->Branch("g_z_prt", &m_z_prt);
  m_trackTree->Branch("px_prt", &m_px_prt);
  m_trackTree->Branch("py_prt", &m_py_prt);
  m_trackTree->Branch("pz_prt", &m_pz_prt);
  m_trackTree->Branch("eta_prt", &m_eta_prt);
  m_trackTree->Branch("pT_prt", &m_pT_prt);
  
  m_trackTree->Branch("nFiltered", &m_nFiltered);
  m_trackTree->Branch("filtered", &m_flt);
  m_trackTree->Branch("eLOC0_flt", &m_eLOC0_flt);
  m_trackTree->Branch("eLOC1_flt", &m_eLOC1_flt);
  m_trackTree->Branch("ePHI_flt", &m_ePHI_flt);
  m_trackTree->Branch("eTHETA_flt", &m_eTHETA_flt);
  m_trackTree->Branch("eQOP_flt", &m_eQOP_flt);
  m_trackTree->Branch("eT_flt", &m_eT_flt);
  m_trackTree->Branch("res_eLOC0_flt", &m_res_eLOC0_flt);
  m_trackTree->Branch("res_eLOC1_flt", &m_res_eLOC1_flt);
  m_trackTree->Branch("res_ePHI_flt", &m_res_ePHI_flt);
  m_trackTree->Branch("res_eTHETA_flt", &m_res_eTHETA_flt);
  m_trackTree->Branch("res_eQOP_flt", &m_res_eQOP_flt);
  m_trackTree->Branch("res_eT_flt", &m_res_eT_flt);
  m_trackTree->Branch("err_eLOC0_flt", &m_err_eLOC0_flt);
  m_trackTree->Branch("err_eLOC1_flt", &m_err_eLOC1_flt);
  m_trackTree->Branch("err_ePHI_flt", &m_err_ePHI_flt);
  m_trackTree->Branch("err_eTHETA_flt", &m_err_eTHETA_flt);
  m_trackTree->Branch("err_eQOP_flt", &m_err_eQOP_flt);
  m_trackTree->Branch("err_eT_flt", &m_err_eT_flt);
  m_trackTree->Branch("pull_eLOC0_flt", &m_pull_eLOC0_flt);
  m_trackTree->Branch("pull_eLOC1_flt", &m_pull_eLOC1_flt);
  m_trackTree->Branch("pull_ePHI_flt", &m_pull_ePHI_flt);
  m_trackTree->Branch("pull_eTHETA_flt", &m_pull_eTHETA_flt);
  m_trackTree->Branch("pull_eQOP_flt", &m_pull_eQOP_flt);
  m_trackTree->Branch("pull_eT_flt", &m_pull_eT_flt);
  m_trackTree->Branch("g_x_flt", &m_x_flt);
  m_trackTree->Branch("g_y_flt", &m_y_flt);
  m_trackTree->Branch("g_z_flt", &m_z_flt);
  m_trackTree->Branch("px_flt", &m_px_flt);
  m_trackTree->Branch("py_flt", &m_py_flt);
  m_trackTree->Branch("pz_flt", &m_pz_flt);
  m_trackTree->Branch("eta_flt", &m_eta_flt);
  m_trackTree->Branch("pT_flt", &m_pT_flt);
  m_trackTree->Branch("chi2", &m_chi2);
  
  m_trackTree->Branch("nSmoothed", &m_nSmoothed);
  m_trackTree->Branch("smoothed", &m_smt);
  m_trackTree->Branch("eLOC0_smt", &m_eLOC0_smt);
  m_trackTree->Branch("eLOC1_smt", &m_eLOC1_smt);
  m_trackTree->Branch("ePHI_smt", &m_ePHI_smt);
  m_trackTree->Branch("eTHETA_smt", &m_eTHETA_smt);
  m_trackTree->Branch("eQOP_smt", &m_eQOP_smt);
  m_trackTree->Branch("eT_smt", &m_eT_smt);
  m_trackTree->Branch("res_eLOC0_smt", &m_res_eLOC0_smt);
  m_trackTree->Branch("res_eLOC1_smt", &m_res_eLOC1_smt);
  m_trackTree->Branch("res_ePHI_smt", &m_res_ePHI_smt);
  m_trackTree->Branch("res_eTHETA_smt", &m_res_eTHETA_smt);
  m_trackTree->Branch("res_eQOP_smt", &m_res_eQOP_smt);
  m_trackTree->Branch("res_eT_smt", &m_res_eT_smt);
  m_trackTree->Branch("err_eLOC0_smt", &m_err_eLOC0_smt);
  m_trackTree->Branch("err_eLOC1_smt", &m_err_eLOC1_smt);
  m_trackTree->Branch("err_ePHI_smt", &m_err_ePHI_smt);
  m_trackTree->Branch("err_eTHETA_smt", &m_err_eTHETA_smt);
  m_trackTree->Branch("err_eQOP_smt", &m_err_eQOP_smt);
  m_trackTree->Branch("err_eT_smt", &m_err_eT_smt);
  m_trackTree->Branch("pull_eLOC0_smt", &m_pull_eLOC0_smt);
  m_trackTree->Branch("pull_eLOC1_smt", &m_pull_eLOC1_smt);
  m_trackTree->Branch("pull_ePHI_smt", &m_pull_ePHI_smt);
  m_trackTree->Branch("pull_eTHETA_smt", &m_pull_eTHETA_smt);
  m_trackTree->Branch("pull_eQOP_smt", &m_pull_eQOP_smt);
  m_trackTree->Branch("pull_eT_smt", &m_pull_eT_smt);
  m_trackTree->Branch("g_x_smt", &m_x_smt);
  m_trackTree->Branch("g_y_smt", &m_y_smt);
  m_trackTree->Branch("g_z_smt", &m_z_smt);
  m_trackTree->Branch("px_smt", &m_px_smt);
  m_trackTree->Branch("py_smt", &m_py_smt);
  m_trackTree->Branch("pz_smt", &m_pz_smt);
  m_trackTree->Branch("eta_smt", &m_eta_smt);
  m_trackTree->Branch("pT_smt", &m_pT_smt);
  */
}
