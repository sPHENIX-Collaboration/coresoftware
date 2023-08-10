#include "CentralityReco.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <calobase/TowerInfoDefs.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <cassert>
#include <vector>
#include <sstream>
#include <string>
#include <TF1.h>
#include <iostream>
#include <fstream>
using namespace std;

CentralityReco::CentralityReco(const std::string& name, const std::string &hist_name, const std::string &tree_name)
  : SubsysReco(name)
  , _verbose(0)
  , _op_mode(0)
  , _mbd_charge_threshold(100.)
  , _zdc_energy_threshold(40.)
{
  _zdc_gain_factors[0] = 1.37;
  _zdc_gain_factors[1] = 0.64;
  _zdc_gain_factors[2] = 0.44;
  _zdc_gain_factors[3] = 1.39;
  _zdc_gain_factors[4] = 0.78;
  _zdc_gain_factors[5] = 0.29;

  _offset = 28.52;
  float mbd_vertex_cuts[5] = {50, 30, 20, 10, 5};
  for (int i = 0; i < 5; i++) _mbd_vertex_cuts[i] = mbd_vertex_cuts[i];

  _hist_filename = hist_name;
  _tree_filename = tree_name;
  std::cout << "CentralityReco::CentralityReco" <<std::endl;  
}

CentralityReco::~CentralityReco()
{
  delete hm;
}

int CentralityReco::Init(PHCompositeNode*)
{
  // Histograms


  fill_n(gaincorr,128,1);
  fill_n(tq_t0_offsets,128,1);

  std::ifstream gainfile ( "/sphenix/user/dlis/Projects/centrality/coresoftware/offline/packages/centrality/gainfile.calib");

  int ch;
  float integ, integerr;
  float peak, peakerr;
  float mbd_width, widtherr;
  float chi2ndf;

  while ( gainfile >> ch >> integ >> peak >> mbd_width >> integerr >> peakerr >> widtherr >> chi2ndf ) {
    gaincorr[ch] = 1.0/peak;
  }
  
  gainfile.close();
  
  std::ifstream tfile ( "/sphenix/user/dlis/Projects/centrality/coresoftware/offline/packages/centrality/bbc_tq_t0.calib");
  
  int pmtnum;
  float meanerr;
  float sigma;
  float sigmaerr;
  for (int ipmt = 0; ipmt < 128; ipmt++) {
    tfile >> pmtnum >> tq_t0_offsets[ipmt] >> meanerr >> sigma >> sigmaerr;

  }

  // ttree
  ttree = new TTree("T","a perservering date tree");

  ttree->Branch("mbd_charge_sum", &_mbd_charge_sum);
  ttree->Branch("mbd_charge_sum_n", &_mbd_charge_sum_n);
  ttree->Branch("mbd_charge_sum_s", &_mbd_charge_sum_s);
  ttree->Branch("mbd_ring_charge_sum_n", _mbd_ring_charge_sum_n, "mbd_ring_charge_sum_n[3]/F");
  ttree->Branch("mbd_ring_charge_sum_s", _mbd_ring_charge_sum_s, "mbd_ring_charge_sum_s[3]/F");

  ttree->Branch("zdc_energy_sum",&_zdc_energy_sum);
  ttree->Branch("zdc_energy_sum_n", &_zdc_energy_sum_n);
  ttree->Branch("zdc_energy_sum_s", &_zdc_energy_sum_s);

  ttree->Branch("mbd_charge", &m_mbd_charge);
  ttree->Branch("mbd_time", &m_mbd_time);
  ttree->Branch("mbd_side", &m_mbd_side);
  ttree->Branch("mbd_channel", &m_mbd_channel);

  ttree->Branch("zdc_energy_low", &m_zdc_energy_low);
  ttree->Branch("zdc_energy_high", &m_zdc_energy_high);
  ttree->Branch("zdc_sum_low", &m_zdc_sum_low);
  ttree->Branch("zdc_sum_high", &m_zdc_sum_high);

  // histograms

  hm = new Fun4AllHistoManager("CENTRALITY_HIST");

  h_mbd_vertex = new TH1D("h_mbd_vertex","",2000, -500, 500); 
  hm->registerHisto(h_mbd_vertex);


  h_mbd_vertex_w_zdc_cut = new TH1D("h_mbd_vertex_zdc_cut","",2000, -500, 500); 
  hm->registerHisto(h_mbd_vertex_w_zdc_cut);

  h_mbd_charge_ns = new TH1D("h_mbd_charge_ns","", 2500, 0, 2500);
  h_mbd_charge_n = new TH1D("h_mbd_charge_n","", 2500, 0, 2500);
  h_mbd_charge_s = new TH1D("h_mbd_charge_s","", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns);
  hm->registerHisto(h_mbd_charge_n);
  hm->registerHisto(h_mbd_charge_s);


  h_mbd_charge_ns_w_zdc_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut","", 2500, 0, 2500);
  h_mbd_charge_n_w_zdc_cut = new TH1D("h_mbd_charge_n_w_zdc_cut","", 2500, 0, 2500);
  h_mbd_charge_s_w_zdc_cut = new TH1D("h_mbd_charge_s_w_zdc_cut","", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns_w_zdc_cut);
  hm->registerHisto(h_mbd_charge_n_w_zdc_cut);
  hm->registerHisto(h_mbd_charge_s_w_zdc_cut);

  h_mbd_charge_ns_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut","", 2500, 0, 2500);
  h_mbd_charge_n_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_n_w_zdc_cut_w_mbd_cut","", 2500, 0, 2500);
  h_mbd_charge_s_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_s_w_zdc_cut_w_mbd_cut","", 2500, 0, 2500);
  hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut);
  for (int i = 0; i < 5; i++)
    {
      h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d",static_cast<int>(_mbd_vertex_cuts[i])),"", 2500, 0, 2500);
      h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex_%d",static_cast<int>(_mbd_vertex_cuts[i])),"", 2500, 0, 2500);
      h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex_%d",static_cast<int>(_mbd_vertex_cuts[i])),"", 2500, 0, 2500);
      hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }
  for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_%d", i),"",2500, 0, 2500);
      h_mbd_ring_charge_sum_s[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_%d", i),"",2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_%d", i),"",2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_%d", i),"",2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_%d", i),"",2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_%d", i),"",2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i]);
      for (int j = 0; j < 5; j++)
	{

	  h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(_mbd_vertex_cuts[j]), i),"",2500, 0, 2500);
	  h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(static_cast<int>(_mbd_vertex_cuts[j])), i),"",2500, 0, 2500);
	  hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);
	  hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);

	}
    }

  // zdc
  h_zdc_energy_ns = new TH1D("h_zdc_energy_ns","",2000, 0, 7000);
  h_zdc_energy_n = new TH1D("h_zdc_energy_n","",2000, 0, 7000);
  h_zdc_energy_s = new TH1D("h_zdc_energy_s","",2000, 0, 7000);
  hm->registerHisto(h_zdc_energy_ns);
  hm->registerHisto(h_zdc_energy_n);
  hm->registerHisto(h_zdc_energy_s);


  // correlations
  h_zdc_mbd_corr_ns = new TH2D("h_zdc_mbd_corr_ns","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n = new TH2D("h_zdc_mbd_corr_n","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s = new TH2D("h_zdc_mbd_corr_s","",100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns);
  hm->registerHisto(h_zdc_mbd_corr_n);
  hm->registerHisto(h_zdc_mbd_corr_s);

  h_zdc_mbd_corr_ns_w_zdc_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n_w_zdc_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s_w_zdc_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut","",100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut);
  hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut);
  hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut);


  h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut","",100, 0, 2500, 100, 0, 10000);
  h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut","",100, 0, 2500, 100, 0, 10000);
  hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut);
  hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut);

  for (int i = 0; i < 5; i++)
    {
      h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])),"",100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])),"",100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])),"",100, 0, 2500, 100, 0, 10000);
      hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }


  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::InitRun(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  CreateNodes(topNode);
  return 0;
}

void CentralityReco::ResetVars()
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  _z_vertex = 0;
  _quality = 0;
  _tubes_hit[0] = 0;
  _tubes_hit[1] = 0;
  _mbd_charge_sum = 0. ;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;
  _zdc_energy_sum = 0.;
  _zdc_energy_sum_n = 0.;
  _zdc_energy_sum_s = 0.;
  for (int i = 0; i < 3; i++)
    {
      _mbd_ring_charge_sum_n[i] = 0.;
      _mbd_ring_charge_sum_s[i] = 0.;
    }
  m_mbd_charge.clear();
  m_mbd_time.clear();
  m_mbd_side.clear();
  m_mbd_channel.clear();

  m_zdc_energy_low.clear();
  m_zdc_sum_low.clear();
  m_zdc_energy_high.clear();
  m_zdc_sum_high.clear();

  return;
}

int CentralityReco::FillVars()
{
  unsigned int size;
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  size = _towers_mbd->size();

  if (size != 256)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  for (unsigned int i = 0; i < size; i++)
    {
      _tmp_tower = _towers_mbd->get_tower_at_channel(i);
      _key = TowerInfoDefs::encode_mbd(i);
      _type = TowerInfoDefs::get_mbd_type(_key);
      if (_type)
	{
	  _energy = _tmp_tower->get_energy();
	  m_mbd_charge.push_back(_energy);
	}
      else
	{
	  _energy = _tmp_tower->get_energy();
	  _side = TowerInfoDefs::get_mbd_side(_key);
	  _channel = TowerInfoDefs::get_mbd_channel(_key);
	  m_mbd_time.push_back(_energy);
	  m_mbd_side.push_back(_side);
	  m_mbd_channel.push_back(_channel);
	  
	}
    }

  size = _towers_zdc->size();
  
  for (unsigned int i = 0; i < size; i++)
    {

      _tmp_tower = _towers_zdc->get_tower_at_channel(i);

      _energy = _tmp_tower->get_energy();

      if (!(i%2))
	{
	  if ((i/2)%4 == 3)
	    {
	      m_zdc_sum_low.push_back(_energy);
	    }
	  else
	    {
	      m_zdc_energy_low.push_back(_energy);
	    }
	}
      else {
	if ((i/2)%4 == 3)
	  {
	    m_zdc_sum_high.push_back(_energy);
	  }
	else
	  {
	    m_zdc_energy_high.push_back(_energy);
	  }
      }
    }
  if (_verbose)
    {
      std::cout << "--------- ZDC data: ----------"<<std::endl;
      std::cout << "South:"<<std::endl;
      for (int i = 0; i < 3; i++)
	{
	  std::cout << i << " : " << m_zdc_energy_low.at(i) << " ("<<m_zdc_energy_high.at(i)<<") "<<std::endl;
	}
	  std::cout << "Sum : " << m_zdc_sum_low.at(0) << " ("<<m_zdc_sum_high.at(0)<<") "<<std::endl;
      std::cout << "North:"<<std::endl;
      for (int i = 0; i < 3; i++)
	{
	  std::cout << i << " : " << m_zdc_energy_low.at(i+3) << " ("<<m_zdc_energy_high.at(i+3)<<") "<<std::endl;
	}
	  std::cout << "Sum : " << m_zdc_sum_low.at(1) << " ("<<m_zdc_sum_high.at(1)<<") "<<std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;

}

void CentralityReco::FillHistograms()
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  h_mbd_vertex->Fill(_z_vertex);

  h_mbd_charge_ns->Fill(_mbd_charge_sum);
  h_mbd_charge_n->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s->Fill(_mbd_charge_sum_s);
  
  for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n[i]->Fill(_mbd_ring_charge_sum_n[i]);
      h_mbd_ring_charge_sum_s[i]->Fill(_mbd_ring_charge_sum_s[i]);
    }
  // zdc
  h_zdc_energy_ns->Fill(_zdc_energy_sum);
  h_zdc_energy_n->Fill(_zdc_energy_sum_n);
  h_zdc_energy_s->Fill(_zdc_energy_sum_s);

  // correlations
  h_zdc_mbd_corr_ns->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);
  
  if (!_zdc_check) return;

  h_mbd_vertex_w_zdc_cut->Fill(_z_vertex);

  h_mbd_charge_ns_w_zdc_cut->Fill(_mbd_charge_sum);
  h_mbd_charge_n_w_zdc_cut->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s_w_zdc_cut->Fill(_mbd_charge_sum_s);
  
  for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n_w_zdc_cut[i]->Fill(_mbd_ring_charge_sum_n[i]);
      h_mbd_ring_charge_sum_s_w_zdc_cut[i]->Fill(_mbd_ring_charge_sum_s[i]);
    }

  // correlations
  h_zdc_mbd_corr_ns_w_zdc_cut->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n_w_zdc_cut->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s_w_zdc_cut->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);

  if (_tubes_hit[0] < 2 || _tubes_hit[1] < 2) return;

  h_mbd_charge_ns_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum);
  h_mbd_charge_n_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_n);
  h_mbd_charge_s_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_s);
  
  for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i]->Fill(_mbd_ring_charge_sum_n[i]);
      h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i]->Fill(_mbd_ring_charge_sum_s[i]);
    }

  // correlations
  h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum, _zdc_energy_sum);
  h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
  h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);

  for (int j = 0; j < 5; j++)
    {
      if (TMath::Abs(_z_vertex) > _mbd_vertex_cuts[j]) continue;
      h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum);
      h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_n);
      h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_s);
  
      for (int i = 0; i < 3; i++)
	{
	  h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j]->Fill(_mbd_ring_charge_sum_n[i]);
	  h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j]->Fill(_mbd_ring_charge_sum_s[i]);
	}

      // correlations
      h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum, _zdc_energy_sum);
      h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_n, _zdc_energy_sum_n);
      h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[j]->Fill(_mbd_charge_sum_s, _zdc_energy_sum_s);
 
    }

  return;
}

void CentralityReco::PrintCentiles()
{

  std::cout << " --------- CentralityReco::PrintCentiles --------- " << std::endl;
  
  float totalevents = h_mbd_charge_ns->Integral();
  int ipercentbin = 1;
  for (int i = 0; i < h_mbd_charge_ns->GetNbinsX();i++)
    {
      if (ipercentbin*.05*totalevents < h_mbd_charge_ns->Integral(0, i + 1))
	{
	  std::cout << " " << ipercentbin*0.05*100 << "% --> " << h_mbd_charge_ns->GetBinCenter(i+1) <<std::endl;
	  ipercentbin++;

	}
    }

  return;
}

int CentralityReco::GetMBDVertexAndCharge()
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  float tdc[2] = {0., 0.};
  unsigned int size = _towers_mbd->size();  
  int side;

  for (unsigned int i = 0; i < size/2; i++)
    {

      side = m_mbd_side.at(i);
      if (m_mbd_charge.at(i) < _mbd_charge_threshold) continue;

      tdc[side] += (25. - m_mbd_time.at(i)*(9.0/5000.) - tq_t0_offsets[m_mbd_channel.at(i) + side*64]);
      
      _tubes_hit[side]++;
					     
    }

  _z_vertex = _offset + 15*(tdc[1]/_tubes_hit[1] - tdc[0]/_tubes_hit[0]);
  
  if (_verbose) std::cout << "Z-vertex = "<<_z_vertex << std::endl;
  for (unsigned int i = 0; i < size/2; i++)
    {
      side = m_mbd_side.at(i);
      
    if (side)
	{
	  _mbd_charge_sum_n += gaincorr[i]*m_mbd_charge.at(i);
	  _mbd_ring_charge_sum_n[mbd_ring_index[i%64]] += gaincorr[i]*m_mbd_charge.at(i);
	}
      else
	{
	  _mbd_charge_sum_s += gaincorr[i]*m_mbd_charge.at(i);
	  _mbd_ring_charge_sum_n[mbd_ring_index[i%64]] += gaincorr[i]*m_mbd_charge.at(i);
	}
    }

  _mbd_charge_sum = _mbd_charge_sum_s + _mbd_charge_sum_n;;
  if (_verbose) std::cout << "mbd charge = "<<_mbd_charge_sum << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::CheckZDC()
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ << std::endl;
  for (unsigned int i = 0; i < 6; i++)
    {
      if (i/3) 
	{
	  _zdc_energy_sum_n += _zdc_gain_factors[i]*m_zdc_energy_low.at(i);
      	  _zdc_energy_sum +=  _zdc_gain_factors[i]*m_zdc_energy_low.at(i);
	}
      else 
	{
	  _zdc_energy_sum_s +=  _zdc_gain_factors[i]*m_zdc_energy_low.at(i);
	  _zdc_energy_sum +=  _zdc_gain_factors[i]*m_zdc_energy_low.at(i);

	}
    }
  if (_verbose) {
    std::cout << " North: "<<_zdc_energy_sum_n << std::endl;
    std::cout << " South: "<<_zdc_energy_sum_s << std::endl;
  }
  _zdc_check = (_zdc_energy_sum_n > _zdc_energy_threshold) && (_zdc_energy_sum_s > _zdc_energy_threshold);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::process_event(PHCompositeNode* topNode)
{
  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  // Get Nodes from the Tree
  if (GetNodes(topNode))
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }


  // Reset Arrays
  ResetVars();

  // Fill Arrays
  if(FillVars())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // Calculate the vertex
  if (GetMBDVertexAndCharge())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // Check the ZDC coincidence
  if (CheckZDC())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  FillHistograms();
  
  ttree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}


int CentralityReco::GetNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;
  if (_op_mode == OperationMode::Performance)
    {
      _central = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");

      if (!_central)
	{
      std::cout << "no centrality node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
	}
    }

  _towers_mbd = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_MBD");
  
  if (!_towers_mbd)
    {
      std::cout << "no mbd towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _towers_zdc = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERS_ZDC");
  
  if (!_towers_zdc)
    {
      std::cout << "no zdc towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityReco::CreateNodes(PHCompositeNode* topNode)
{

  if (_verbose) std::cout << __FILE__ << " :: "<<__FUNCTION__ <<" :: " << __LINE__ << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    }

  PHNodeIterator dstIter(dstNode);

  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!detNode)
    {
      std::cout << PHWHERE << "Detector Node missing, making one"<<std::endl;
      detNode = new PHCompositeNode("GLOBAL");
      dstNode->addNode(detNode);
    }
  

  if (_op_mode == OperationMode::Performance)
    {
      CentralityInfov1 *central = new CentralityInfov1();
      
      PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, "CentralityInfo", "PHObject");
      detNode->addNode(centralityNode);
    }
  return;
}


int CentralityReco::End(PHCompositeNode* /* topNode*/)
{
  outfile = new TFile(_tree_filename.c_str(), "RECREATE");
  ttree->Write();
  outfile->Close();

  PrintCentiles();
  hm->dumpHistos(_hist_filename.c_str(), "RECREATE");
  return 0;
}
