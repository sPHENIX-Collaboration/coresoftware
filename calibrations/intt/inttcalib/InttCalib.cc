#include "InttCalib.h"

#include <cdbobjects/CDBTTree.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>

#include <boost/format.hpp>

#include <cmath>
#include <iostream>
#include <filesystem>
#include <limits>

InttCalib::InttCalib(const std::string& name)
  : SubsysReco(name)
{
}

int InttCalib::Init()
{
  m_evts = 0;
  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    for (int bco = 0; bco < 129; ++bco)
    {
      m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] = 0;
    }
  }

  m_do_nothing = false;

  if(m_survey.LoadFromCDB("InttSurveyMap"))
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not load 'InttSurveyMap' from CDB\n"
              << "\tModule will do nothing" << std::endl;
	m_do_nothing = true;
    // gSystem->Exit(1);
    // exit(1);
  }

  if (m_feemap.LoadFromCDB("InttFeeMap"))
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not load 'InttFeeMap' from CDB\n"
              << "\tModule will do nothing" << std::endl;
	m_do_nothing = true;
    // gSystem->Exit(1);
    // exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::InitRun(PHCompositeNode* /*unused*/)
{
  return Init();
}

int InttCalib::process_event(PHCompositeNode* top_node)
{
  if(m_do_nothing)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  InttRawHitContainer* intt_raw_hit_container = findNode::getClass<InttRawHitContainer>(top_node, "INTTRAWHIT");
  if (!intt_raw_hit_container)
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not get 'INTTRAWHIT' from node tree\n"
              << "\tModule will do nothing" << std::endl;
	m_do_nothing = true;
    // gSystem->Exit(1);
    // exit(1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  for (size_t n = 0, N = intt_raw_hit_container->get_nhits(); n < N; ++n)
  {
    InttRawHit* intt_raw_hit = intt_raw_hit_container->get_hit(n);
    if (!intt_raw_hit)
    {
		continue;
    }

    InttMap::RawData_s raw{
        .pid = intt_raw_hit->get_packetid(),             //
        .fee = intt_raw_hit->get_fee(),                  //
        .chp = (intt_raw_hit->get_chip_id() + 25) % 26,  //
        .chn = intt_raw_hit->get_channel_id(),           //
    };

    int bco_diff = ((intt_raw_hit->get_bco() & 0x7fU) - intt_raw_hit->get_FPHX_BCO() + 128) % 128;

    ++m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][bco_diff];
    ++m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128];
  }

  ++m_evts;
  if ((Verbosity() > 1) && ((m_evts % 1000) == 0))
  {
    std::cout << "Finished event: " << m_evts << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::EndRun(int const run_number)
{
  if(m_do_nothing)
  {
    std::cout << PHWHERE << "\n"
              << "\tMember 'm_do_nothing' set\n"
              << "\tDoing nothing" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_run_num = run_number;

  ConfigureHotMap_v2();
  MakeHotMapCdb_v2();
  MakeHotMapPng_v2();

  ConfigureBcoMap();
  MakeBcoMapCdb();
  MakeBcoMapPng();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::SaveHitrates()
{
  TFile* file = TFile::Open("/sphenix/user/jbertaux/hitrates.root", "RECREATE");
  if(!file)
  {
    std::cerr << "\n" << PHWHERE << "\n\tfile\n" << std::endl;
	return 1;
  }
  file->cd();
  TTree* tree = new TTree("hitrate_tree", "hitrate_tree");
  tree->SetDirectory(file);

  InttMap::RawData_s raw;
  tree->Branch("pid", &raw.pid);
  tree->Branch("fee", &raw.fee);
  tree->Branch("chp", &raw.chp);
  tree->Branch("chn", &raw.chn);

  double hitrate;
  tree->Branch("hitrate", &hitrate);

  for(raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
	tree->Fill();
  }

  tree->Write();
  file->Write();
  file->Close();

  return 0;
}

int InttCalib::LoadHitrates()
{
  TFile* file = TFile::Open("/sphenix/user/jbertaux/hitrates.root", "READ");
  if(!file)
  {
    std::cerr << "\n" << PHWHERE << "\n\tfile\n" << std::endl;
	return 1;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get("hitrate_tree"));
  if(!tree)
  {
    std::cerr << "\n" << PHWHERE << "\n\ttree\n" << std::endl;
	return 1;
  }

  InttMap::RawData_s raw;
  tree->SetBranchAddress("pid", &raw.pid);
  tree->SetBranchAddress("fee", &raw.fee);
  tree->SetBranchAddress("chp", &raw.chp);
  tree->SetBranchAddress("chn", &raw.chn);

  double hitrate;
  tree->SetBranchAddress("hitrate", &hitrate);

  for(Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n)
  {
    tree->GetEntry(n);
    m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] = hitrate;
  }

  m_evts = 1.0;

  return 0;
}

int InttCalib::ConfigureHotMap_v2()
{
  std::map<double, int> hitrate_pdf[5];
  std::string name[5] = {"inner_a", "inner_b", "outer_a", "outer_b", "entire"};

  for(InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
    InttMap::Offline_s ofl = m_feemap.ToOffline(raw);

    int index = 0;
	index += (ofl.layer < 5) ? 0 : 2; // +2 for outer
	index += (ofl.ladder_z % 2) ? 0 : 1; // +1 for type B

	++hitrate_pdf[index][hitrate];
	++hitrate_pdf[4][hitrate];
  }

  for(int i = 4; i >= 0; --i)
  {
    ConfigureHist (
      m_hist[i],
      m_fit[i],
      hitrate_pdf[i],
      name[i]
    );

	double mean = m_fit[i]->GetParameter(1);
	double sigma = m_fit[i]->GetParameter(2);

	m_min[i] = mean - m_NUM_SIGMA * sigma;
	m_max[i] = mean + m_NUM_SIGMA * sigma;

	if(m_min[i] < 0)m_min[i] = 0;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapCdb_v2()
{
  if (m_hotmap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree* cdbttree = new CDBTTree(m_hotmap_cdb_file);
  int size = 0;
  for(InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
    InttMap::Offline_s ofl = m_feemap.ToOffline(raw);

    int index = 0;
	index += (ofl.layer < 5) ? 0 : 2; // +2 for outer
	index += (ofl.ladder_z % 2) ? 0 : 1; // +1 for type B

	if(m_min[index] < hitrate && hitrate < m_max[index])
	{
      continue;
	}

    cdbttree->SetIntValue(size, "felix_server",  raw.pid - 3001);
    cdbttree->SetIntValue(size, "felix_channel", raw.fee);
    cdbttree->SetIntValue(size, "chip",          raw.chp);
    cdbttree->SetIntValue(size, "channel",       raw.chn);
    ++size;
  }
  cdbttree->SetSingleIntValue("size", size);

  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapPng_v2()
{
  std::string name[5] = {"inner_a", "inner_b", "outer_a", "outer_b", "entire"};

  // make png
  TCanvas* cnvs = new TCanvas (
    "cnvs",   //
    "cnvs",   //
    1280, 720 //
  );

  cnvs->cd();
  TPad* hist_pad = new TPad (
    "hist_pad", "hist_pad",
    0.0, 0.0,
	0.8, 1.0
  );
  hist_pad->SetFillStyle(4000);
  hist_pad->Range(0.0, 0.0, 1.0, 1.0);
  hist_pad->Draw();

  hist_pad->cd();
  for(int i = 4; i >= 0; --i)
  {
	m_hist[i]->SetLineColor(5 - i);
	m_hist[i]->SetLineWidth(2);
	m_hist[i]->Draw("same");

	m_fit[i]->SetLineColor(5 - i);
	m_fit[i]->SetLineWidth(2);
	m_fit[i]->Draw("same");

    double y_max = m_hist[4]->GetBinContent(m_hist[4]->GetMaximumBin());

	TLine line;
	line.SetLineColor(5 - i);
	line.SetLineWidth(2);
    line.DrawLine(m_min[i], 0, m_min[i], y_max);
    line.DrawLine(m_max[i], 0, m_max[i], y_max);
  }

  cnvs->cd();
  TPad* lgnd_pad = new TPad (
    "lgnd_pad", "lgnd_pad",
    0.8, 0.0,
	1.0, 1.0
  );
  lgnd_pad->SetFillStyle(4000);
  lgnd_pad->Range(0.0, 0.0, 1.0, 1.0);
  lgnd_pad->Draw();

  lgnd_pad->cd();
  for(int i = 4; i >= 0; --i)
  {
    TText text;
	text.SetTextColor(5 - i);
	text.SetTextAlign(22);
	text.SetTextSize(0.2);
	text.DrawText(0.5, 0.1 + i * 0.2, name[i].c_str());
  }

  cnvs->Update();
  cnvs->Show();

  return 0;
}

int InttCalib::ConfigureHotMap()
{
  m_hitrates.clear();
  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
    InttMap::Offline_s ofl;
    if (m_feemap.Convert(ofl, raw))
    {
      continue;
    }

    if (adjust_hitrate(ofl, hitrate))
    {
      continue;
    }

    ++m_hitrates[hitrate];
  }

  m_invcdf.clear();
  double total = 0;
  for (auto const& [hitrate, count] : m_hitrates)
  {
    total += count;
    double fraction = total / m_NUM_CHANNELS;
    m_invcdf[fraction] = hitrate;
  }

  double prev_hitrate = 0.0, prev_fraction = 0.0;
  double best_min_ratio = 0.0, best_max_ratio = 0.0;
  for (auto const& [fraction, hitrate] : m_invcdf)
  {
    double ratio = (hitrate - prev_hitrate) / (fraction - prev_fraction);
    prev_hitrate = hitrate;
    prev_fraction = fraction;

    if (best_min_ratio < ratio && fraction < 0.5)
    {
      best_min_ratio = ratio;
      m_min_hitrate = hitrate;
      m_min_fraction = fraction;
    }

    if (best_max_ratio < ratio && 0.5 < fraction)
    {
      best_max_ratio = ratio;
      m_max_hitrate = hitrate;
      m_max_fraction = 1.0 - prev_fraction;
    }

    if (best_min_ratio < best_max_ratio)
    {
      break;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int InttCalib::MakeHotMapCdb()
{
  if (m_hotmap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree* cdbttree = new CDBTTree(m_hotmap_cdb_file);

  // Dummy calibration for now
  // Mask exactly one channel, which we plan on masking at Felix level anyway

  cdbttree->SetIntValue(0, "felix_server", 3002);
  cdbttree->SetIntValue(0, "felix_channel", 2);
  cdbttree->SetIntValue(0, "chip", 14);
  cdbttree->SetIntValue(0, "channel", 0);

  cdbttree->SetSingleIntValue("size", 1);

  // int size = 0;
  // for(InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  // {
  //   double hitrate = (double)m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / (double)m_evts;
  //   InttMap::Offline_s ofl;
  //   if(m_feemap.Convert(ofl, raw))
  //   {
  //      continue;
  //   }

  //   if(adjust_hitrate(ofl, hitrate))
  //   {
  //     continue;
  //   }

  //   if(m_min_hitrate < hitrate && hitrate < m_max_hitrate)
  //   {
  //      continue;
  //   }

  //   cdbttree->SetIntValue(size, "felix_server",  raw.pid - 3001);
  //   cdbttree->SetIntValue(size, "felix_channel", raw.fee);
  //   cdbttree->SetIntValue(size, "chip",          raw.chp);
  //   cdbttree->SetIntValue(size, "channel",       raw.chn);
  //   ++size;
  // }
  // cdbttree->SetSingleIntValue("size", size);

  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapPng()
{
  if (m_hotmap_png_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Canvas
  gStyle->SetOptStat(0);
  TCanvas* hitrate_cnvs = new TCanvas(  //
      "hitrate_cnvs", "hitrate_cnvs",   //
      1280, 720                         //
  );
  hitrate_cnvs->Draw();

  // PDF
  double q1 = 0.0, q2 = 0.0, q3 = 0.0;
  for (auto const& [fraction, hitrate] : m_invcdf)
  {
    if (fraction < 0.25)
    {
      q1 = hitrate;
    }
    if (fraction < 0.5)
    {
      q2 = hitrate;
    }
    if (fraction < 0.75)
    {
      q3 = hitrate;
    }
  }

  // Freedman-Diaconis rule
  double bin_width = 2 * (q3 - q1) / pow(m_NUM_CHANNELS, 1.0 / 3.0);
  double lower = 0.0, upper = q2 + 2.0 * (q3 - q1);
  int num_bins = std::floor((upper - lower) / bin_width);

  hitrate_cnvs->cd();
  TPad* hitrate_pdf_pad = new TPad(          //
      "hitrate_pdf_pad", "hitrate_pdf_pad",  //
      0.0, 0.2, 0.5, 1.0                     //
  );
  hitrate_pdf_pad->SetFillStyle(4000);  // transparent
  hitrate_pdf_pad->Range(0.0, 0.0, 1.0, 1.0);
  hitrate_pdf_pad->Draw();

  hitrate_pdf_pad->cd();
  TH1D* hitrate_pdf_hist = new TH1D(           //
      "hitrate_pdf_hist", "hitrate_pdf_hist",  //
      num_bins, lower, upper                   //
  );
  hitrate_pdf_hist->SetTitle("Hitrate PDF;Adjusted Hitrate;Count");
  hitrate_pdf_hist->Draw();

  // Fill
  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
    InttMap::Offline_s ofl;
    if (m_feemap.Convert(ofl, raw))
    {
      continue;
    }

    if (adjust_hitrate(ofl, hitrate))
    {
      continue;
    }

    hitrate_pdf_hist->Fill(hitrate);
  }
  double y_max = 1.1 * hitrate_pdf_hist->GetBinContent(hitrate_pdf_hist->GetMaximumBin());
  hitrate_pdf_hist->GetYaxis()->SetRangeUser(0, y_max);
  TLine line;
  line.SetLineWidth(1);
  line.SetLineColor(kRed);
  line.DrawLine(m_min_hitrate, 0, m_min_hitrate, y_max);
  line.DrawLine(m_max_hitrate, 0, m_max_hitrate, y_max);

  // CDF
  hitrate_cnvs->cd();
  TPad* hitrate_cdf_pad = new TPad(          //
      "hitrate_cdf_pad", "hitrate_cdf_pad",  //
      0.5, 0.2, 1.0, 1.0                     //
  );
  hitrate_cdf_pad->SetLogy();
  hitrate_cdf_pad->SetFillStyle(4000);  // transparent
  hitrate_cdf_pad->Range(0.0, 0.0, 1.0, 1.0);
  hitrate_cdf_pad->Draw();

  hitrate_cdf_pad->cd();
  TH1D* hitrate_cdf_hist = new TH1D(           //
      "hitrate_cdf_hist", "hitrate_cdf_hist",  //
      500, 0.0, 100.0                          //
  );
  hitrate_cdf_hist->SetTitle("Hitrate Inverse CDF;Percentile;Adjusted Hitrate");
  hitrate_cdf_hist->GetXaxis()->SetNdivisions(10, kTRUE);
  hitrate_cdf_hist->GetYaxis()->SetRangeUser(0.1, 10 * m_hitrates.rbegin()->first);
  hitrate_cdf_hist->Draw();

  for (auto& [fraction, hitrate] : m_invcdf)
  {
    hitrate_cdf_hist->SetBinContent(hitrate_cdf_hist->GetXaxis()->FindBin(100 * fraction), hitrate);
  }
  for (int i = 2; i < hitrate_cdf_hist->GetNbinsX(); ++i)
  {
    if (hitrate_cdf_hist->GetBinContent(i))
    {
      continue;
    }
    hitrate_cdf_hist->SetBinContent(i, hitrate_cdf_hist->GetBinContent(i - 1));
  }
  line.SetLineWidth(1);
  line.SetLineColor(kRed);
  line.DrawLine(0, m_min_hitrate, 100, m_min_hitrate);
  line.DrawLine(0, m_max_hitrate, 100, m_max_hitrate);

  // Caption
  hitrate_cnvs->cd();
  TPad* caption_pad = new TPad(
      "caption_pad", "caption_pad",  //
      0.0, 0.0, 1.0, 0.2             //
  );
  caption_pad->SetFillStyle(4000);  // transparent
  caption_pad->Range(0.0, 0.0, 1.0, 1.0);
  caption_pad->Draw();

  caption_pad->cd();
  TText caption_text;
  caption_text.SetTextAlign(22);
  caption_text.SetTextSize(0.20);
  caption_text.SetTextColor(kBlack);
  caption_text.DrawText(0.5, 0.75, (boost::format("Run: %08d Events: %d") % m_run_num % m_evts).str().c_str());

  caption_text.SetTextSize(0.10);
  caption_text.DrawText(0.5, 0.5, (boost::format("%.3E <= [hitrate] excludes %06.3lf%%") % m_min_hitrate % (m_min_fraction * 100)).str().c_str());
  caption_text.DrawText(0.5, 0.35, (boost::format("[hitrate] <= %.3E excludes %06.3lf%%") % m_max_hitrate % (m_max_fraction * 100)).str().c_str());
  caption_text.DrawText(0.5, 0.2, (boost::format("Keeping %.3lf%%") % (100.0 - (m_min_fraction + m_max_fraction) * 100.0)).str().c_str());

  hitrate_cnvs->Update();
  hitrate_cnvs->Show();
  hitrate_cnvs->SaveAs(m_hotmap_png_file.c_str());

  delete hitrate_pdf_hist;
  delete hitrate_cdf_hist;
  delete hitrate_cnvs;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::ConfigureBcoMap()
{
  m_bcorates.clear();
  for (InttMap::RawData_s raw{.pid = InttMap::RawDataBegin.pid, .fee = InttMap::RawDataBegin.fee}; raw != InttMap::RawDataEnd; ++raw)
  {
    for (auto& bco_count : m_bcorates[raw])
    {
      bco_count = 0;
    }
  }

  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / m_evts;
    InttMap::Offline_s ofl;
    if (m_feemap.Convert(ofl, raw))
    {
      continue;
    }

    if (adjust_hitrate(ofl, hitrate))
    {
      continue;
    }

    if (hitrate <= m_min_hitrate || m_max_hitrate <= hitrate)
    {
      continue;
    }

    for (int bco = 0; bco < 128; ++bco)
    {
      m_bcorates[raw][bco] += m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][bco];
    }
  }

  m_bcopeaks.clear();
  for (auto const& [raw, bco_arr] : m_bcorates)
  {
    int max_counts = 0, bco_peak = 0;
    for (int bco = 0; bco < 128; ++bco)
    {
      if (max_counts < bco_arr[bco])
      {
        bco_peak = bco;
        max_counts = bco_arr[bco];
      }
    }
    m_bcopeaks[raw] = bco_peak;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeBcoMapCdb()
{
  if (m_bcomap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree* cdbttree = new CDBTTree(m_bcomap_cdb_file);

  int size = 0;
  for (auto const& [raw, bco] : m_bcopeaks)
  {
    cdbttree->SetIntValue(size, "felix_server", raw.pid - 3001);
    cdbttree->SetIntValue(size, "felix_channel", raw.fee);
    cdbttree->SetIntValue(size, "bco_diff", bco);
    ++size;
  }
  cdbttree->SetSingleIntValue("size", size);

  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeBcoMapPng()
{
  if (m_bcomap_png_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Canvas
  gStyle->SetOptStat(0);
  TCanvas* bco_cnvs = new TCanvas(  //
      "bco_cnvs", "bco_cnvs",       //
      1280, 720                     //
  );
  bco_cnvs->Draw();

  TH1D* bco_hist[112] = {};
  for (int i = 0; i < 16; ++i)
  {
    std::string tpadtitle = (boost::format("bco_pad_%02d") % i).str();
    bco_cnvs->cd();
    TPad* bco_pdf_pad = new TPad(  //
        tpadtitle.c_str(), tpadtitle.c_str(),
        // NOLINTNEXTLINE(bugprone-integer-division)
        (i % 4 + 0.0) / 4.0 * 0.9 + 0.0, (3.0 - i / 4) / 4.0 * 0.9 + 0.1,
        // NOLINTNEXTLINE(bugprone-integer-division)
        (i % 4 + 1.0) / 4.0 * 0.9 + 0.0, (4.0 - i / 4) / 4.0 * 0.9 + 0.1);
    bco_pdf_pad->SetFillStyle(4000);  // transparent
    bco_pdf_pad->SetLogy();
    bco_pdf_pad->SetLeftMargin(0.15);
    bco_pdf_pad->SetRightMargin(0.05);
    if ((i / 4) % 2)
    {
      bco_pdf_pad->SetTopMargin(0.0);
      bco_pdf_pad->SetBottomMargin(0.15);
    }
    else
    {
      bco_pdf_pad->SetTopMargin(0.15);
      bco_pdf_pad->SetBottomMargin(0.0);
    }
    bco_pdf_pad->Range(0.0, 0.0, 1.0, 1.0);
    bco_pdf_pad->Draw();

    InttMap::RawData_s raw_begin{.pid = (i % 4) + 4 * (i / 8) + 3001, .fee = (i / 4) % 2 ? 0 : 7};
    InttMap::RawData_s raw_end{.pid = (i % 4) + 4 * (i / 8) + 3001, .fee = (i / 4) % 2 ? 6 : 13};

    double max = 0;
    for (InttMap::RawData_s raw = raw_begin; raw <= raw_end; ++raw)
    {
      int h = (raw.pid - 3001) * 14 + raw.fee;
      bco_pdf_pad->cd();
      std::string htitle = (boost::format("bco_hist_%03d") % h).str();
      bco_hist[h] = new TH1D(              //
          htitle.c_str(), htitle.c_str(),  //
          128, -0.5, 127.5                 //
      );

      bco_hist[h]->SetTitle((boost::format(";BCO Difference;intt%01d (%01d - %02d)") % (raw_begin.pid - 3001) % raw_begin.fee % raw_end.fee).str().c_str());
      bco_hist[h]->GetYaxis()->CenterTitle();
      bco_hist[h]->GetYaxis()->SetTitleSize(0.12);
      bco_hist[h]->GetYaxis()->SetTitleOffset(0.6);
      bco_hist[h]->GetYaxis()->SetLabelSize(0.07);
      bco_hist[h]->GetXaxis()->SetTitleSize(0.10);
      bco_hist[h]->GetXaxis()->SetTitleOffset(0.6);
      bco_hist[h]->GetXaxis()->SetLabelSize(0.07);
      bco_hist[h]->SetLineColor(GetFeeColor(raw.fee));
      bco_hist[h]->Draw("same");

      for (int bco = 0; bco < 128; ++bco)
      {
        bco_hist[h]->SetBinContent(bco + 1, m_bcorates[raw][bco]);
        if (max < m_bcorates[raw][bco])
        {
          max = m_bcorates[raw][bco];
        }
      }
    }
    for (InttMap::RawData_s raw = raw_begin; raw <= raw_end; ++raw)
    {
      int h = (raw.pid - 3001) * 14 + raw.fee;
      bco_hist[h]->GetYaxis()->SetRangeUser(0.1, 10 * max);
    }
  }

  // Legend
  bco_cnvs->cd();
  TPad* legend_pad = new TPad(
      "legend_pad", "legend_pad",  //
      0.9, 0.1, 1.0, 1.0           //
  );
  legend_pad->SetFillStyle(4000);  // transparent
  legend_pad->Range(0.0, 0.0, 1.0, 1.0);
  legend_pad->Draw();
  legend_pad->cd();

  for (int fee = 0; fee < 7; ++fee)
  {
    TText legend_text;
    legend_text.SetTextAlign(22);
    legend_text.SetTextColor(kBlack);
    legend_text.SetTextSize(0.1);
    legend_text.DrawText(0.6, (2.0 * fee + 1.0) / 14.0, (boost::format("FCh %01d, %02d") % fee % (fee + 7)).str().c_str());

    double x[4] = {-1.0, +1.0, +1.0, -1.0};
    double y[4] = {-1.0, -1.0, +1.0, +1.0};
    for (int i = 0; i < 4; ++i)
    {
      x[i] *= 0.1;
      x[i] += 0.2;

      y[i] *= 0.008;
      y[i] += (2.0 * fee + 1.0) / 14.0;
    }

    TPolyLine box;
    box.SetFillColor(GetFeeColor(fee));
    box.SetLineColor(kBlack);
    box.SetLineWidth(1);
    box.DrawPolyLine(4, x, y, "f");
  }

  // Caption
  bco_cnvs->cd();
  TPad* caption_pad = new TPad(
      "caption_pad", "caption_pad",  //
      0.0, 0.0, 1.0, 0.1             //
  );
  caption_pad->SetFillStyle(4000);  // transparent
  caption_pad->Range(0.0, 0.0, 1.0, 1.0);
  caption_pad->Draw();

  caption_pad->cd();
  TText caption_text;
  caption_text.SetTextAlign(22);
  caption_text.SetTextSize(0.20);
  caption_text.SetTextColor(kBlack);
  caption_text.DrawText(0.5, 0.75, (boost::format("Run: %08d Events: %d") % m_run_num % m_evts).str().c_str());

  bco_cnvs->Update();
  bco_cnvs->Show();
  bco_cnvs->SaveAs(m_bcomap_png_file.c_str());

  for (auto& hist : bco_hist)
  {
    delete hist;
  }
  delete bco_cnvs;

  return Fun4AllReturnCodes::EVENT_OK;
}

void InttCalib::Debug()
{
  Init();
  LoadHitrates();

  ConfigureHotMap_v2();
  // MakeHotMapCdb_v2();
  MakeHotMapPng_v2();
}

int InttCalib::ConfigureHist(TH1D*& hist, TF1*& fit, std::map<double, int> const& hitrate_map, std::string const& name)
{
  // quartiles (less sentive to outliers, better to configure with)
  double N_entries = 0.0;
  double quartile[5] = {};
  for(auto const& [hitrate, count] : hitrate_map)
  {
    N_entries += count;
  }

  double sum = 0.0;
  for(auto const& [hitrate, count] : hitrate_map)
  {
    for(int i = 0; i < 5; ++i)
    {
      if(sum / N_entries < 0.25 * i)
      {
        quartile[i] = hitrate;
      }
    }
    sum += count;
  }

  double lower = 0;
  double upper = quartile[2] + 3.0 * (quartile[3] - quartile[1]);
  int n_edges = 2;
  for(auto const& [hitrate, count] : hitrate_map)
  {
    if(hitrate <= lower)continue;
    if(upper <= hitrate)continue;
    ++n_edges;
  }
  double* bins = new double[n_edges];
  bins[0] = lower;
  bins[n_edges - 1] = upper;
  n_edges = 1;
  for(auto const& [hitrate, count] : hitrate_map)
  {
    if(hitrate <= lower)continue;
    if(upper <= hitrate)continue;
    bins[n_edges] = hitrate;
	++n_edges;
  }

  // Freedman-Diaconis rule
  //   https://en.wikipedia.org/wiki/Freedman-Diaconis_rule
  // double bin_width = 2.0 * (quartile[3] - quartile[1]) / pow(N_entries, 1.0 / 3.0);
  // int N_bins = std::floor((upper - lower) / bin_width) + 1;

  if(Verbosity())
  {
    std::cout << "size: " << hitrate_map.size() << std::endl;
    std::cout << "N_entries: " << N_entries << std::endl;
	std::cout << "quartiles: " << std::endl;
    for(auto const& q : quartile)
    {
      std::cout << "\t" << q << std::endl;
    }
  }

  // Make hist
  delete hist;
  hist = new TH1D (           //
    (name + "_hist").c_str(), //
    (name + "_hist").c_str(), //
	n_edges - 1, bins         //
  );
  delete[] bins;

  for(auto const& [hitrate, count] : hitrate_map)
  {
    for(int i = 0; i < count; ++i)
	{
      hist->Fill(hitrate);
	}
  }

  delete fit;
  fit = new TF1 (            //
    (name + "_fit").c_str(), //
	"gaus",                  //
	lower, upper             //
  );

  fit->SetParameter(0, N_entries);                   // normalization
  fit->SetParameter(1, quartile[2]);                 // mean ~ median
  fit->SetParameter(2, (quartile[3] - quartile[1])); // standard deviation ~ IQR

  if(Verbosity())
  {
    hist->Fit(fit, "RNL"); // range, no-draw, log likelihood
  }
  else
  {
    hist->Fit(fit, "RNLQ"); // range, no-draw, log likelihood, quiet
  }

  return 0;
}

int InttCalib::adjust_hitrate(InttMap::Offline_s const& ofl, double& hitrate) const
{
  hitrate /= (ofl.ladder_z % 2) ? 2.0 : 1.6;  // normalize by sensor length
  InttSurveyMap::val_t transform;
  if (m_survey.GetStripTransform(ofl, transform))
  {
    std::cout << PHWHERE << "\n"
              << "InttSurveyMap::GetStripTransform failed for\n"
              << ofl
              << std::endl;
    return 1;
  }
  Eigen::Vector3d v = transform.translation() - m_vertex;
  Eigen::Vector3d n{transform.linear()(0, 1), transform.linear()(1, 1), transform.linear()(2, 1)};
  hitrate *= v.squaredNorm();
  hitrate /= (-1.0 * n.dot(v) / v.norm());
  return 0;
}

Color_t InttCalib::GetFeeColor(int fee) const
{
  switch (fee % 7)
  {
  case 1:
    return kRed;
  case 2:
    return kGreen;
  case 3:
    return kYellow;
  case 4:
    return kBlue;
  case 5:
    return kMagenta;
  case 6:
    return kCyan;
  default:  // 0
    break;
  }
  return kBlack;
}
