#include "CentralityReco.h"

#include "CentralityInfov2.h"
#include "bbc/BbcDefs.h"
#include "bbc/BbcVertexMapv1.h"
#include "bbc/BbcVertexv2.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

CentralityReco::CentralityReco(const std::string &name, const std::string &hist_name, const std::string &tree_name)
  : SubsysReco(name)
{
  _zdc_gain_factors[0] = 1.37;
  _zdc_gain_factors[1] = 0.64;
  _zdc_gain_factors[2] = 0.44;
  _zdc_gain_factors[3] = 1.39;
  _zdc_gain_factors[4] = 0.78;
  _zdc_gain_factors[5] = 0.29;

  const int centrality_map[20] = {1999, 1499, 1291, 1102, 937, 790, 660, 547, 449, 363, 289, 227, 174, 130, 94, 66, 45, 0, 0, 0};
  for (int i = 0; i < 20; i++)
  {
    _centrality_map[i] = centrality_map[i];
  }

  _offset = 28.52;

  const float mbd_vertex_cuts[5] = {50, 30, 20, 10, 5};
  for (int i = 0; i < 5; i++)
  {
    _mbd_vertex_cuts[i] = mbd_vertex_cuts[i];
  }

  _hist_filename = hist_name;
  _tree_filename = tree_name;
}

CentralityReco::~CentralityReco()
{
  delete hm;
}

int CentralityReco::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  std::fill_n(gaincorr, 128, 1);
  std::fill_n(tq_t0_offsets, 128, 1);
  const char *gaincalib = getenv("CENTRALITY_GAINCALIB");
  std::string gainfilename;
  if (gaincalib == nullptr)
  {
    const char *offline_main = getenv("OFFLINE_MAIN");
    assert(offline_main);  // make cppcheck happy
    gainfilename = offline_main;
    gainfilename += "/share/centrality/gainfile.calib";
  }
  else
  {
    gainfilename = gaincalib;
  }
  if (!std::filesystem::exists(gainfilename))
  {
    std::cout << PHWHERE << gainfilename << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  std::ifstream gainfile(gainfilename);

  int ch;
  float integ;
  float integerr;
  float peak;
  float peakerr;
  float mbd_width;
  float widtherr;
  float chi2ndf;

  while (gainfile >> ch >> integ >> peak >> mbd_width >> integerr >> peakerr >> widtherr >> chi2ndf)
  {
    gaincorr[ch] = 1.0 / peak;
  }

  gainfile.close();

  const char *bbc_tq_t0calib = getenv("CENTRALITY_BBC_TQ_T0CALIB");
  std::string bbc_tq_t0_filename;
  if (bbc_tq_t0calib == nullptr)
  {
    const char *offline_main = getenv("OFFLINE_MAIN");
    assert(offline_main);  // make cppcheck happy
    bbc_tq_t0_filename = offline_main;
    bbc_tq_t0_filename += "/share/centrality/bbc_tq_t0.calib";
  }
  else
  {
    bbc_tq_t0_filename = bbc_tq_t0calib;
  }
  if (!std::filesystem::exists(bbc_tq_t0_filename))
  {
    std::cout << PHWHERE << bbc_tq_t0_filename << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  std::ifstream tfile(bbc_tq_t0_filename);

  int pmtnum;
  float meanerr;
  float sigma;
  float sigmaerr;
  for (float &tq_t0_offset : tq_t0_offsets)
  {
    tfile >> pmtnum >> tq_t0_offset >> meanerr >> sigma >> sigmaerr;
  }
  if (_op_mode == OperationMode::QA)
  {
    outfile = new TFile(_tree_filename.c_str(), "RECREATE");
    // ttree
    ttree = new TTree("T", "a perservering date tree");

    ttree->Branch("z_vertex", &_z_vertex);
    ttree->Branch("isMinBias", &_isMinBias);
    ttree->Branch("mbd_charge_sum", &_mbd_charge_sum);
    ttree->Branch("mbd_charge_sum_n", &_mbd_charge_sum_n);
    ttree->Branch("mbd_charge_sum_s", &_mbd_charge_sum_s);
    ttree->Branch("mbd_ring_charge_sum_n", _mbd_ring_charge_sum_n, "mbd_ring_charge_sum_n[3]/F");
    ttree->Branch("mbd_ring_charge_sum_s", _mbd_ring_charge_sum_s, "mbd_ring_charge_sum_s[3]/F");

    ttree->Branch("zdc_energy_sum", &_zdc_energy_sum);
    ttree->Branch("zdc_energy_sum_n", &_zdc_energy_sum_n);
    ttree->Branch("zdc_energy_sum_s", &_zdc_energy_sum_s);

    ttree->Branch("mbd_charge", m_mbd_charge, "mbd_charge[128]/F");
    ttree->Branch("mbd_time", m_mbd_time, "mbd_time[128]/F");
    ttree->Branch("mbd_charge_raw", m_mbd_charge_raw, "mbd_charge_raw[128]/F");
    ttree->Branch("mbd_time_raw", m_mbd_time_raw, "mbd_time_raw[128]/F");

    ttree->Branch("mbd_side", m_mbd_side, "mbd_side[128]/I");
    ttree->Branch("mbd_channel", m_mbd_channel, "mbd_channel[128]/I");

    ttree->Branch("zdc_energy_low", m_zdc_energy_low, "zdc_energy_low[6]/F");
    ttree->Branch("zdc_energy_high", m_zdc_energy_high, "zdc_energy_high[6]/F");
    ttree->Branch("zdc_sum_low", m_zdc_sum_low, "zdc_sum_low[2]/F");
    ttree->Branch("zdc_sum_high", m_zdc_sum_high, "zdc_sum_high[2]/F");

    // histograms

    hm = new Fun4AllHistoManager("CENTRALITY_HIST");

    h_mbd_vertex = new TH1D("h_mbd_vertex", "", 2000, -500, 500);
    hm->registerHisto(h_mbd_vertex);

    h_mbd_vertex_w_zdc_cut = new TH1D("h_mbd_vertex_zdc_cut", "", 2000, -500, 500);
    hm->registerHisto(h_mbd_vertex_w_zdc_cut);

    h_mbd_charge_ns = new TH1D("h_mbd_charge_ns", "", 2500, 0, 2500);
    h_mbd_charge_n = new TH1D("h_mbd_charge_n", "", 2500, 0, 2500);
    h_mbd_charge_s = new TH1D("h_mbd_charge_s", "", 2500, 0, 2500);
    hm->registerHisto(h_mbd_charge_ns);
    hm->registerHisto(h_mbd_charge_n);
    hm->registerHisto(h_mbd_charge_s);

    h_mbd_charge_ns_w_zdc_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut", "", 2500, 0, 2500);
    h_mbd_charge_n_w_zdc_cut = new TH1D("h_mbd_charge_n_w_zdc_cut", "", 2500, 0, 2500);
    h_mbd_charge_s_w_zdc_cut = new TH1D("h_mbd_charge_s_w_zdc_cut", "", 2500, 0, 2500);
    hm->registerHisto(h_mbd_charge_ns_w_zdc_cut);
    hm->registerHisto(h_mbd_charge_n_w_zdc_cut);
    hm->registerHisto(h_mbd_charge_s_w_zdc_cut);

    h_mbd_charge_ns_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
    h_mbd_charge_n_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_n_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
    h_mbd_charge_s_w_zdc_cut_w_mbd_cut = new TH1D("h_mbd_charge_s_w_zdc_cut_w_mbd_cut", "", 2500, 0, 2500);
    hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut);
    hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut);
    hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut);
    for (int i = 0; i < 5; i++)
    {
      h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH1D(Form("h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }
    for (int i = 0; i < 3; i++)
    {
      h_mbd_ring_charge_sum_n[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut[i]);

      h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_%d", i), "", 2500, 0, 2500);
      h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_%d", i), "", 2500, 0, 2500);
      hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[i]);
      hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[i]);
      for (int j = 0; j < 5; j++)
      {
        h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(_mbd_vertex_cuts[j]), i), "", 2500, 0, 2500);
        h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j] = new TH1D(Form("h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex_%d_%d", static_cast<int>(static_cast<int>(_mbd_vertex_cuts[j])), i), "", 2500, 0, 2500);
        hm->registerHisto(h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);
        hm->registerHisto(h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[i][j]);
      }
    }

    // zdc
    h_zdc_energy_ns = new TH1D("h_zdc_energy_ns", "", 2000, 0, 7000);
    h_zdc_energy_n = new TH1D("h_zdc_energy_n", "", 2000, 0, 7000);
    h_zdc_energy_s = new TH1D("h_zdc_energy_s", "", 2000, 0, 7000);
    hm->registerHisto(h_zdc_energy_ns);
    hm->registerHisto(h_zdc_energy_n);
    hm->registerHisto(h_zdc_energy_s);

    // correlations
    h_zdc_mbd_corr_ns = new TH2D("h_zdc_mbd_corr_ns", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_n = new TH2D("h_zdc_mbd_corr_n", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_s = new TH2D("h_zdc_mbd_corr_s", "", 100, 0, 2500, 100, 0, 10000);
    hm->registerHisto(h_zdc_mbd_corr_ns);
    hm->registerHisto(h_zdc_mbd_corr_n);
    hm->registerHisto(h_zdc_mbd_corr_s);

    h_zdc_mbd_corr_ns_w_zdc_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_n_w_zdc_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_s_w_zdc_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut", "", 100, 0, 2500, 100, 0, 10000);
    hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut);
    hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut);
    hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut);

    h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
    h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut = new TH2D("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut", "", 100, 0, 2500, 100, 0, 10000);
    hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut);
    hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut);
    hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut);

    for (int i = 0; i < 5; i++)
    {
      h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i] = new TH2D(Form("h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex_%d", static_cast<int>(_mbd_vertex_cuts[i])), "", 100, 0, 2500, 100, 0, 10000);
      hm->registerHisto(h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[i]);
      hm->registerHisto(h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[i]);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  CreateNodes(topNode);
  return 0;
}

void CentralityReco::ResetVars()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _z_vertex = 0;
  _quality = 0;
  _isMinBias = 0;
  _tubes_hit[0] = 0;
  _tubes_hit[1] = 0;
  _mbd_charge_sum = 0.;
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
  for (int i = 0; i < 128; i++)
  {
    m_mbd_charge[i] = 0.;
    m_mbd_time[i] = 0.;
    m_mbd_charge_raw[i] = 0.;
    m_mbd_time_raw[i] = 0.;
    m_mbd_side[i] = 0;
    m_mbd_channel[i] = 0;
  }
  for (int i = 0; i < 6; i++)
  {
    m_zdc_energy_low[i] = 0;
    m_zdc_energy_high[i] = 0;
  }
  for (int i = 0; i < 2; i++)
  {
    m_zdc_sum_low[i] = 0;
    m_zdc_sum_high[i] = 0;
  }
  return;
}

int CentralityReco::FillVars()
{
  unsigned int size;
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
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
      _side = TowerInfoDefs::get_mbd_side(_key);
      _channel = TowerInfoDefs::get_mbd_channel(_key);
      _energy = _tmp_tower->get_energy();
      if (_type)
	{
	  m_mbd_charge_raw[_channel + _side * 64] = _energy;
	  m_mbd_charge[_channel + _side * 64] = gaincorr[_channel + _side * 64] * _energy;
	}
      else
	{
	  m_mbd_time_raw[_channel + _side * 64] = _energy;
	  m_mbd_time[_channel + _side * 64] = (25. - _energy * (9.0 / 5000.) - tq_t0_offsets[_channel + _side * 64]);
	  m_mbd_side[_channel + _side * 64] = _side;
	  m_mbd_channel[_channel + _side * 64] = _channel;
	}
    }

  size = _towers_zdc->size();

  for (unsigned int i = 0; i < size; i++)
  {
    _tmp_tower = _towers_zdc->get_tower_at_channel(i);
    _energy = _tmp_tower->get_energy();

    if (!(i % 2))
    {
      if ((i / 2) % 4 == 3)
      {
        m_zdc_sum_low[i / 8] = _energy;
      }
      else
      {
        m_zdc_energy_low[i / 2] = _zdc_gain_factors[i / 2] * _energy;
      }
    }
    else
    {
      if ((i / 2) % 4 == 3)
      {
        m_zdc_sum_high[i / 8] = _energy;
      }
      else
      {
        m_zdc_energy_high[i / 2] = _zdc_gain_factors[i / 2] * _energy;
      }
    }
  }

  if (Verbosity() > 5)
  {
    std::cout << "--------- ZDC data: ----------" << std::endl;
    std::cout << "South:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
      std::cout << i << " : " << m_zdc_energy_low[i] << " (" << m_zdc_energy_high[i] << ") " << std::endl;
    }
    std::cout << "Sum : " << m_zdc_sum_low[0] << " (" << m_zdc_sum_high[0] << ") " << std::endl;
    std::cout << "North:" << std::endl;
    for (int i = 0; i < 3; i++)
    {
      std::cout << i << " : " << m_zdc_energy_low[i + 3] << " (" << m_zdc_energy_high[i + 3] << ") " << std::endl;
    }
    std::cout << "Sum : " << m_zdc_sum_low[1] << " (" << m_zdc_sum_high[1] << ") " << std::endl;
  }
  if (Verbosity() > 5)
  {
    std::cout << "--------- MBD data: ----------" << std::endl;
    std::cout << "South:" << std::endl;
    for (int i = 0; i < 64; i++)
    {
      std::cout << m_mbd_channel[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time[i] << ") " << std::endl;
    }
    std::cout << "North:" << std::endl;
    for (int i = 64; i < 128; i++)
    {
      std::cout << m_mbd_channel[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time[i] << ") " << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityReco::FillHistograms()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
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

  if (!_zdc_check)
  {
    return;
  }

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
  if (Verbosity())
  {
    std::cout << "tubes hit in N/S : " << _tubes_hit[1] << " / " << _tubes_hit[0] << std::endl;
  }
  if (_tubes_hit[0] < 2 || _tubes_hit[1] < 2)
  {
    return;
  }

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
    if (std::abs(_z_vertex) > _mbd_vertex_cuts[j])
    {
      continue;
    }
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
  for (int i = 0; i < h_mbd_charge_ns->GetNbinsX(); i++)
  {
    if (ipercentbin * .05 * totalevents < h_mbd_charge_ns->Integral(0, i + 1))
    {
      std::cout << " " << ipercentbin * 0.05 * 100 << "% --> " << h_mbd_charge_ns->GetBinCenter(i + 1) << std::endl;
      ipercentbin++;
    }
  }

  return;
}

int CentralityReco::GetMBDVertexAndCharge()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _tdc[0] = 0.;
  _tdc[1] = 0.;
  unsigned int size = _towers_mbd->size();
  int side;

  for (unsigned int i = 0; i < size / 2; i++)
  {
    side = m_mbd_side[i];
    if (m_mbd_charge[i] <= _mbd_charge_threshold)
    {
      continue;
    }

    _tdc[side] += m_mbd_time[i];  //(25. - m_mbd_time[i]*(9.0/5000.) - tq_t0_offsets[m_mbd_channel[i] + side*64]);

    _tubes_hit[side]++;
  }

  if ((_tubes_hit[0] >= 2) && (_tubes_hit[1] >= 2))
  {
    _isMinBias = 1;
  }
  else
  {
    _isMinBias = 0;
  }

  _z_vertex = _offset + 15 * (_tdc[1] / _tubes_hit[1] - _tdc[0] / _tubes_hit[0]);

  if (Verbosity())
  {
    std::cout << "Z-vertex = " << _z_vertex << std::endl;
  }
  for (unsigned int i = 0; i < size / 2; i++)
  {
    side = m_mbd_side[i];

    if (side)
    {
      _mbd_charge_sum_n += m_mbd_charge[i];
      _mbd_ring_charge_sum_n[mbd_ring_index[i % 64]] += m_mbd_charge[i];
    }
    else
    {
      _mbd_charge_sum_s += m_mbd_charge[i];
      _mbd_ring_charge_sum_s[mbd_ring_index[i % 64]] += m_mbd_charge[i];
    }
  }

  _mbd_charge_sum = _mbd_charge_sum_s + _mbd_charge_sum_n;
  
  if (Verbosity())
  {
    std::cout << "mbd charge = " << _mbd_charge_sum << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::CheckZDC()
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  for (unsigned int i = 0; i < 6; i++)
  {
    if (i / 3)
    {
      _zdc_energy_sum_n += m_zdc_energy_low[i];
      _zdc_energy_sum += m_zdc_energy_low[i];
    }
    else
    {
      _zdc_energy_sum_s += m_zdc_energy_low[i];
      _zdc_energy_sum += m_zdc_energy_low[i];
    }
  }
  if (Verbosity())
  {
    std::cout << " North: " << _zdc_energy_sum_n << std::endl;
    std::cout << " South: " << _zdc_energy_sum_s << std::endl;
  }
  _zdc_check = (_zdc_energy_sum_n > _zdc_energy_threshold) && (_zdc_energy_sum_s > _zdc_energy_threshold);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias
  bool final_check = _zdc_check && _isMinBias;
  _central->setMinBias(final_check);

  // Fill z-vertex

  auto vtx = std::make_unique<BbcVertexv2>();

  vtx->set_z(_z_vertex);
  vtx->set_t((_tdc[0] + _tdc[1])/2.);
  vtx->set_bbc_ns(0, _tubes_hit[0], _mbd_charge_sum_s, _tdc[0]);
  vtx->set_bbc_ns(1, _tubes_hit[1], _mbd_charge_sum_n, _tdc[1]);

  _bbc_vertex_map->insert(vtx.release());
  
  // Fill Centrality

  float value = 0;
  for (int i = 0; i < 20; i++)
  {
    if (_centrality_map[i] < _mbd_charge_sum)
    {
      value = 0.05 * i;
      break;
    }
  }
  if (!value)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _central->set_centile(CentralityInfo::PROP::mbd_NS, value);

  if (Verbosity()) _bbc_vertex_map->identify();
  if (Verbosity()) _bbc_vertex_map->begin()->second->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Reset Arrays
  ResetVars();

  // Fill Arrays
  if (FillVars())
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

  if (_op_mode == OperationMode::Performance)
  {
    if (FillCentralityInfo())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if (Verbosity())
    {
      _central->identify();
    }
  }
  else if (_op_mode == OperationMode::QA)
  {
    FillHistograms();
    ttree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  _bbc_vertex_map = findNode::getClass<BbcVertexMapv1>(topNode, "BbcVertexMap");
  
  if (!_bbc_vertex_map)
    {
      std::cout << "no mbd vertex node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  if (_op_mode == OperationMode::Performance)
  {
    _central = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");

    if (!_central)
    {
      std::cout << "no centrality node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  _towers_mbd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_MBD");
  
  if (!_towers_mbd)
    {
      std::cout << "no mbd towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");

  if (!_towers_zdc)
  {
    std::cout << "no zdc towers node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityReco::CreateNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

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
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

  if (_op_mode == OperationMode::Performance)
  {
    CentralityInfov2 *central = new CentralityInfov2();

    PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, "CentralityInfo", "PHObject");
    detNode->addNode(centralityNode);
  }

  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "BBC"));
  if (!bbcNode)
  {
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    bbcNode = new PHCompositeNode("BBC");
    dstNode->addNode(bbcNode);
  }

  _bbc_vertex_map = findNode::getClass<BbcVertexMapv1>(topNode, "BbcVertexMap");
  if (!_bbc_vertex_map)
    {
      _bbc_vertex_map = new BbcVertexMapv1();
      PHIODataNode<PHObject> *vertexNode = new PHIODataNode<PHObject>(_bbc_vertex_map, "BbcVertexMap", "PHObject");
      detNode->addNode(vertexNode);
    }

  return;
}

int CentralityReco::End(PHCompositeNode * /* topNode*/)
{
  if (_op_mode == OperationMode::QA)
  {
    outfile->Write();
    outfile->Close();

    PrintCentiles();
    hm->dumpHistos(_hist_filename.c_str(), "RECREATE");
  }
  return 0;
}
