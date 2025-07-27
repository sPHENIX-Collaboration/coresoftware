#include "InttCalib.h"

#include <intt/InttMapping.h>

#include <cdbobjects/CDBTTree.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>  // for size_t
#include <cstdint>  // for uint64_t
#include <sstream>  // for basic_ostringstream

#include <format>
#include <iostream>
#include <limits>

InttCalib::InttCalib(const std::string &name)
  : SubsysReco(name)
{
}

InttCalib::~InttCalib()
{
  for (auto &hist : m_hist)
  {
    delete hist;
  }
  for (auto &fit : m_fit)
  {
    delete fit;
  }
  for (auto &hist : m_hist_fee)
  {
    delete hist;
  }
  for (auto &fit : m_fit_fee)
  {
    delete fit;
  }
  for (auto &hist : m_hist_half)
  {
    delete hist;
  }
  for (auto &hist : m_bco_peak)
  {
    delete hist;
  }
}

int InttCalib::InitRun(PHCompositeNode * /*unused*/)
{
  m_evts = 0;
  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    for (int bco = 0; bco < 129; ++bco)
    {
      m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][bco] = 0;
    }
  }

  m_do_nothing = false;

  std::cout << "INITRUNEND" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::process_event(PHCompositeNode *top_node)
{
  if (m_do_nothing)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(top_node, "GL1RAWHIT");
  uint64_t gl1_bco = (gl1 != nullptr) ? (gl1->getBCO() & 0xFFFFFFFFFFU)
                                      : std::numeric_limits<uint64_t>::max();

  if (!gl1)
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not get 'GL1RAWHIT' from node tree" << std::endl;
    if (m_streaming)
    {
      std::cout << "\tRunmode is Streaming \n"
                << "\tModule will do nothing" << std::endl;
      m_do_nothing = true;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    std::cout << "\tRunmode is Triggered \n"
              << "\tModule will process without GL1" << std::endl;
  }
  InttRawHitContainer *intt_raw_hit_container =
      findNode::getClass<InttRawHitContainer>(top_node, m_rawhit_container_name);
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
    InttRawHit *intt_raw_hit = intt_raw_hit_container->get_hit(n);
    if (!intt_raw_hit)
    {
      continue;
    }

    InttNameSpace::RawData_s raw{
        .felix_server = intt_raw_hit->get_packetid() - 3001,  //
        .felix_channel = intt_raw_hit->get_fee(),             //
        .chip = (intt_raw_hit->get_chip_id() + 25) % 26,      //
        .channel = intt_raw_hit->get_channel_id(),            //
    };

    // Trigger mode BCO Offset
    // old convention
    // int bco_diff = ((intt_raw_hit->get_bco() & 0x7fU) -
    // intt_raw_hit->get_FPHX_BCO() + 128) % 128; int bco_diff =
    // (intt_raw_hit->get_FPHX_BCO() - (intt_raw_hit->get_bco() & 0x7fU) + 128)
    // % 128;
    int bco_diff =
        (m_streaming) ?
                      // Streaming mode BCO Offset : Hit BCO(FPHX BCO + INTT
                      // BCO) - GL1 BCO
            intt_raw_hit->get_FPHX_BCO() + intt_raw_hit->get_bco() - gl1_bco
                      :
                      // Trigger mode BCO Offset
            (intt_raw_hit->get_FPHX_BCO() - (intt_raw_hit->get_bco() & 0x7fU) +
             128) %
                128;
    // Not use the hit from abort gap region
    if (m_streaming && ((intt_raw_hit->get_FPHX_BCO() > 116) ||
                        (intt_raw_hit->get_FPHX_BCO() < 5)))
    {
      bco_diff = -999;
    }

    if (bco_diff > -1)
    {
      ++m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][bco_diff];
    }
    ++m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128];
  }

  ++m_evts;
  if ((Verbosity() > 1) && ((m_evts % 1000) == 0))
  {
    std::cout << "event: " << m_evts << std::endl;
  }
  if (m_evts == m_evts_bco && m_evts_bco != 0)
  {
    ConfigureBcoMap();
    MakeBcoMapCdb();
    MakeBcoMapPng();
    m_do_make_bco = false;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::EndRun(int const run_number)
{
  if (m_do_nothing)
  {
    std::cout << PHWHERE << "\n"
              << "\tMember 'm_do_nothing' set\n"
              << "\tDoing nothing" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_run_num = run_number;
  if (m_do_fee)
  {
    ConfigureHotMap_fee();
    MakeHotMapCdb_fee();
    MakeHotMapROOT_fee();
  }
  else
  {
    ConfigureHotMap_v3();
    MakeHotMapCdb_v3();
    MakeHotMapPng_v3();
  }
  if (m_do_make_bco)
  {
    ConfigureBcoMap();
    MakeBcoMapCdb();
    MakeBcoMapPng();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int InttCalib::ConfigureHotMap_fee()
{
  std::map<double, int> hitrate_pdf[m_MAX_LADDER]{};
  std::string name[m_MAX_LADDER];
  std::string title[m_MAX_LADDER];
  for (int i = 0; i < m_MAX_LADDER; ++i)
  {
    // name[i] = std::format("intt{:01d}", (i / 4));
    name[i] = std::format("h_InttCalib_intt{:01d}_fee{}", (i / 14), (i % 14));
    title[i] = name[i];
  }

  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    double hitrate =
        m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

    int index = GetFeeIndex(raw, ofl);
    adjust_hitrate(ofl, hitrate);

    ++hitrate_pdf[index][hitrate];
  }
  std::vector<double> middle_keys(m_MAX_LADDER);
  for (int i = 0; i < m_MAX_LADDER; ++i)
  {
    size_t map_size = hitrate_pdf[i].size();

    size_t mid_index = map_size / 2;
    double middle_key = 0.;

    auto it = hitrate_pdf[i].begin();
    std::advance(it, mid_index);
    middle_key = it->first;
    middle_keys[i] = middle_key;
  }
  std::sort(middle_keys.begin(), middle_keys.end());
  double global_maxbin = 5 * (middle_keys[56] + middle_keys[57]) / 2.0;
  for (int i = 0; i < m_MAX_LADDER; ++i)
  {
    ConfigureHist_v3(m_hist_fee[i], m_fit_fee[i], global_maxbin, hitrate_pdf[i], name[i], title[i]);

    double mean = m_fit_fee[i]->GetParameter(1);
    double sigma = m_fit_fee[i]->GetParameter(2);
    m_mean_fee[i] = mean;
    m_sigma_fee[i] = sigma;

    m_min_fee[i] = mean - m_NUM_SIGMA_COLD * sigma;
    if (m_min_fee[i] <= 0)
    {
      m_min_fee[i] = -999;
    }
    m_max_fee[i] = mean + m_NUM_SIGMA_HOT * sigma;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::ConfigureHotMap_v3()
{
  std::map<double, int> hitrate_pdf[m_MAX_INDEX]{};
  std::string name[m_MAX_INDEX];
  std::string title[m_MAX_INDEX];
  for (int i = 0; i < m_MAX_INDEX; ++i)
  {
    // name[i] = std::format("intt{:01d}", % (i / 4));
    name[i] = std::format("h_InttCalib_intt{:01d}", i);
    title[i] = name[i];
  }

  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != raw.felix_server)
    {
      continue;
    }
    double hitrate =
        m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

    int index = GetIndex(raw, ofl);
    adjust_hitrate(ofl, hitrate);

    ++hitrate_pdf[index][hitrate];
  }
  std::vector<double> middle_keys(8);
  double global_maxbin = 0;
  for (int i = 0; i < m_MAX_INDEX; ++i)
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != i)
    {
      continue;
    }
    size_t map_size = hitrate_pdf[i].size();

    size_t mid_index = map_size / 2;
    double middle_key = 0.;

    auto it = hitrate_pdf[i].begin();
    std::advance(it, mid_index);
    middle_key = it->first;
    middle_keys[i] = middle_key;
    global_maxbin = 5 * middle_key;
  }
  std::sort(middle_keys.begin(), middle_keys.end());
  if (m_FELIX_TARGET == -1)
  {
    global_maxbin = 5 * (middle_keys[3] + middle_keys[4]) / 2.0;
  }
  for (int i = 0; i < m_MAX_INDEX; ++i)
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != i)
    {
      continue;
    }
    ConfigureHist_v3(m_hist[i], m_fit[i], global_maxbin, hitrate_pdf[i], name[i], title[i]);
    QAHistManagerDef::getHistoManager()->registerHisto(m_hist[i]);
    QAHistManagerDef::getHistoManager()->registerHisto(m_fit[i]);
    int nBins = m_hist[i]->GetNbinsX();
    double xMin = m_hist[i]->GetXaxis()->GetXmin();
    double xMax = m_hist[i]->GetXaxis()->GetXmax();

    m_hist_half[i] =
        new TH1D(std::format("h_InttCalib_half_hist_{}", i).c_str(),
                 "New Histogram with Same Binning", nBins, xMin, xMax);

    double mean = m_fit[i]->GetParameter(1);
    double sigma = m_fit[i]->GetParameter(2);
    m_mean[i] = mean;
    m_sigma[i] = sigma;
    // m_min[i] = mean - m_NUM_SIGMA * sigma;

    m_min[i] = mean - m_NUM_SIGMA_COLD * sigma;
    if (m_min[i] <= 0)
    {
      m_min[i] = -999;
    }
    m_max[i] = mean + m_NUM_SIGMA_HOT * sigma;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapCdb_fee()
{
  if (m_hotmap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree *cdbttree = new CDBTTree(m_hotmap_cdb_file);
  int size = 0;
  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    double hitrate =
        m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

    int index = GetFeeIndex(raw, ofl);
    adjust_hitrate(ofl, hitrate);
    // Example Part How to add channel masking mannually.
    //  if( (raw.pid-3001 == 2) && (raw.felix_channel == 9) && (raw.chip == 15))
    //  {
    //    cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
    //    cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
    //    cdbttree->SetIntValue(size, "chip", raw.chip);
    //    cdbttree->SetIntValue(size, "channel", raw.channel);
    //    cdbttree->SetIntValue(size, "flag", 4);
    //    ++size;
    //    continue;
    //  }
    // End of Example Part
    if (hitrate == 0)
    {  // dead channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 1);
      ++size;
      continue;
    }
    if (m_min_fee[index] < hitrate && hitrate < m_max_fee[index])
    {  // good channel
      continue;
    }
    if (hitrate > m_max_fee[index])
    {  // hot channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 8);
      ++size;
      continue;
    }
    if (hitrate < m_min_fee[index])
    {  // cold channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 4);
      ++size;
      continue;
    }
  }
  cdbttree->SetSingleIntValue("size", size);
  cdbttree->SetSingleIntValue("event", m_evts);
  for (int i = 0; i < m_MAX_LADDER; i++)
  {
    std::string meanname = "mean" + std::to_string(i);
    std::string sigmaname = "sigma" + std::to_string(i);
    cdbttree->SetSingleDoubleValue(meanname, m_mean_fee[i]);
    cdbttree->SetSingleDoubleValue(sigmaname, m_sigma_fee[i]);
  }

  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapCdb_v3()
{
  if (m_hotmap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree *cdbttree = new CDBTTree(m_hotmap_cdb_file);
  int size = 0;
  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != raw.felix_server)
    {
      continue;
    }
    double hitrate =
        m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

    int index = GetIndex(raw, ofl);
    adjust_hitrate(ofl, hitrate);
    if (m_half_min[index] < hitrate && hitrate < m_half_max[index])
    {
      m_hitmap_half[raw.felix_server][raw.felix_channel][raw.chip]++;
      // if (m_hitmap_half[raw.felix_server][raw.felix_channel][raw.chip] > 100)
      // {
      //   cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      //   cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      //   cdbttree->SetIntValue(size, "chip", raw.chip);
      //   cdbttree->SetIntValue(size, "channel", raw.channel);
      //   cdbttree->SetIntValue(size, "flag", 2);
      //   ++size;
      //   continue;
      // }
    }

    if (hitrate == 0)
    {  // dead channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 1);
      ++size;
      continue;
    }
    if (m_min[index] < hitrate && hitrate < m_max[index])
    {  // good channel
      continue;
    }
    if (hitrate > m_max[index])
    {  // hot channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 8);
      ++size;
      continue;
    }
    if (hitrate < m_min[index])
    {  // cold channel
      cdbttree->SetIntValue(size, "felix_server", raw.felix_server);
      cdbttree->SetIntValue(size, "felix_channel", raw.felix_channel);
      cdbttree->SetIntValue(size, "chip", raw.chip);
      cdbttree->SetIntValue(size, "channel", raw.channel);
      cdbttree->SetIntValue(size, "flag", 4);
      ++size;
      continue;
    }
  }
  cdbttree->SetSingleIntValue("size", size);
  cdbttree->SetSingleIntValue("event", m_evts);
  for (int i = 0; i < 8; i++)
  {
    std::string meanname = "mean" + std::to_string(i);
    std::string sigmaname = "sigma" + std::to_string(i);
    cdbttree->SetSingleDoubleValue(meanname, m_mean[i]);
    cdbttree->SetSingleDoubleValue(sigmaname, m_sigma[i]);
  }

  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapROOT_fee()
{
  const int rows = 8;
  const int cols = 14;
  const int numHists = rows * cols;
  TCanvas *c1 = new TCanvas("c1", "Histograms with Fits", 2000, 1400);
  c1->Divide(cols, rows);

  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.01);
  gStyle->SetPadLeftMargin(0.01);
  gStyle->SetPadRightMargin(0.01);

  for (int i = 0; i < numHists; ++i)
  {
    c1->cd(i + 1);
    m_hist_fee[i]->Draw();
    m_fit_fee[i]->Draw("same");
  }
  if (!m_hotmap_png_file.empty())
  {
    c1->SaveAs(m_hotmap_png_file.c_str());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeHotMapPng_v3()
{
  // canvas
  gStyle->SetOptStat(0);
  TCanvas *cnvs = new TCanvas("hitrate_cnvs",  //
                              "hitrate_cnvs",  //
                              1280, 720        //
  );
  double n_hot = 0;
  double n_cold = 0;
  double n_dead = 0;
  double n_total = 0;
  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    double hitrate =
        m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    InttNameSpace::Offline_s ofl = InttNameSpace::ToOffline(raw);

    int index = GetIndex(raw, ofl);
    adjust_hitrate(ofl, hitrate);

    if (!(m_min[index] < hitrate))
    {
      ++n_cold;
    }
    if (!(hitrate < m_max[index]))
    {
      ++n_hot;
    }
    if (hitrate == 0)
    {
      ++n_dead;
    }
    ++n_total;
    if (m_hitmap_half[raw.felix_server][raw.felix_channel][raw.chip] > 100)
    {
      m_hist_half[raw.felix_server]->Fill(hitrate);
    }
  }
  for (int j = 0; j < 8; ++j)
  {
    std::string name = std::format("hist_pad_{:01d}", j);

    cnvs->cd();
    TPad *hist_pad = new TPad(  //
        name.c_str(),
        name.c_str(),  //
                       // (j % 4 + 0.0) / 4.0 * 0.9 + 0.0, (1.0 - j / 4) / 2.0 *
                       // 0.9 + 0.1, // (j % 4 + 1.0) / 4.0 * 0.9 + 0.0, (2.0 - j
                       // / 4) / 2.0 * 0.9 + 0.1  //
                       // NOLINTNEXTLINE(bugprone-integer-division)
        ((j % 4 + 0.0) / 4.0 * 1.0) + 0.0,
        // NOLINTNEXTLINE(bugprone-integer-division)
        ((1.0 - j / 4) / 2.0 * 0.9) + 0.1,  //
                                            // NOLINTNEXTLINE(bugprone-integer-division)
        ((j % 4 + 1.0) / 4.0 * 1.0) + 0.0,
        // NOLINTNEXTLINE(bugprone-integer-division)
        ((2.0 - j / 4) / 2.0 * 0.9) + 0.1  //
    );

    hist_pad->SetFillStyle(4000);
    hist_pad->Range(0.0, 0.0, 1.0, 1.0);
    hist_pad->SetLogy();
    hist_pad->Draw();

    hist_pad->cd();
    double x_max = 0;
    double y_max = 0;
    // for(int i = j * 4; i < (j + 1) * 4; ++i)
    for (int i = j; i < j + 1; ++i)
    {
      // m_hist[i]->SetLineColor(GetFeeColor(i - j * 4));
      if (!m_hist[i])
      {
        continue;
      }
      m_hist[i]->SetLineColor(kBlack);
      m_hist[i]->SetLineWidth(2);
      m_hist_half[i]->SetLineColor(kRed);
      m_hist_half[i]->SetLineWidth(3);

      // m_fit[i]->SetLineColor(GetFeeColor(i - j * 4));
      m_fit[i]->SetLineColor(kBlue);
      m_fit[i]->SetLineWidth(2);

      double temp_max;

      temp_max = m_hist[i]->GetBinContent(m_hist[i]->GetMaximumBin());
      y_max = std::max(y_max, temp_max);

      temp_max = m_hist[i]->GetXaxis()->GetBinLowEdge(
          m_hist[i]->GetXaxis()->GetNbins() - 1);
      temp_max += m_hist[i]->GetXaxis()->GetBinWidth(
          m_hist[i]->GetXaxis()->GetNbins() - 1);
      x_max = std::max(x_max, temp_max);
    }
    y_max *= 10;

    for (int i = j; i < j + 1; ++i)
    {
      if (!m_hist[i])
      {
        continue;
      }
      m_hist[i]->GetXaxis()->SetRangeUser(0, x_max);
      m_hist[i]->GetYaxis()->SetRangeUser(1, y_max);
      m_hist[i]->Draw("same");
      // histogram for half entry chips. Not used for now
      // m_hist_half[i]->Draw("same");
      m_fit[i]->Draw("same");

      TLine line;
      line.SetLineColor(kRed);
      line.SetLineWidth(2);
      line.DrawLine(m_max[i], 0, m_max[i], y_max);

      TLine line2;
      line2.SetLineColor(kBlue);
      line2.DrawLine(m_min[i], 0, m_min[i], y_max);
    }
  }
  cnvs->cd();
  TPad *legend_pad = new TPad("legend_pad", "legend_pad", 0.9, 0.1, 1.0, 1.0);
  legend_pad->SetFillStyle(4000);
  legend_pad->Range(0.0, 0.0, 1.0, 1.0);

  legend_pad->cd();
  for (int i = 0; i < 4; ++i)
  {
    TText text;
    text.SetTextColor(GetFeeColor(i));
    text.SetTextAlign(22);
    text.SetTextSize(0.15);
    if (!m_hist[i])
    {
      continue;
    }
    std::string title = m_hist[i]->GetName();
    text.DrawText(0.5, (2.0 * i + 1.0) / (2.0 * 4), title.substr(6, 7).c_str());
  }

  cnvs->cd();
  TPad *caption_pad =
      new TPad("caption_pad", "caption_pad", 0.0, 0.0, 1.0, 0.1);
  caption_pad->SetFillStyle(4000);
  caption_pad->Range(0.0, 0.0, 1.0, 1.0);
  caption_pad->Draw();

  caption_pad->cd();
  TText caption;
  caption.SetTextColor(kBlack);
  caption.SetTextAlign(22);
  caption.SetTextSize(0.25);
  caption.DrawText(0.5, 0.75, std::format("Run: {:08d} Events: {}", m_run_num, m_evts).c_str());
  caption.DrawText(0.5, 0.50,
                   std::format("Fraction Cold: {:.3f}% Fraction Dead: {:.3f}%", (n_cold * 100 / n_total), (n_dead * 100 / n_total)).c_str());
  caption.DrawText(0.5, 0.25, std::format("Fraction Hot: {:.3f}%", (n_hot * 100 / n_total)).c_str());

  cnvs->Update();
  cnvs->Show();
  if (!m_hotmap_png_file.empty())
  {
    cnvs->SaveAs(m_hotmap_png_file.c_str());
  }

  delete cnvs;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::ConfigureBcoMap()
{
  for (int felix = 0; felix < 8; felix++)
  {
    for (int fee = 0; fee < 14; fee++)
    {
      if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != felix)
      {
        continue;
      }
      // Find chip with highest total hits for this pid/fee across all channels
      int chp_most = 0;
      int max_hits = 0;
      for (int chp = 0; chp < 26; chp++)
      {
        int chip_total = 0;
        // Sum hits across all channels for this chip
        for (int chan = 0; chan < 128; chan++)
        {
          chip_total += m_hitmap[felix][fee][chp][chan][128];
        }
        if (chip_total > max_hits)
        {
          max_hits = chip_total;
          chp_most = chp;
        }
      }
      for (int chp = 0; chp < 26; chp++)
      {
        // Sum hits across all channels for this chip
        for (int chan = 0; chan < 128; chan++)
        {
          if (chp != chp_most)
          {
            for (int bco = 0; bco < 128; ++bco)
            {
              m_bcorates_fee[felix][fee][bco] +=
                  m_hitmap[felix][fee][chp][chan][bco];
            }
          }
        }
      }
    }
  }

  for (int felix = 0; felix < 8; felix++)
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != felix)
    {
      continue;
    }
    for (int fee = 0; fee < 14; fee++)
    {
      int max_counts = 0;
      int bco_peak = 0;
      for (int bco = 0; bco < 128; bco++)
      {
        if (max_counts < m_bcorates_fee[felix][fee][bco])
        {
          bco_peak = bco;
          max_counts = m_bcorates_fee[felix][fee][bco];
        }
      }

      if (max_counts < 50)  // if max_count is less than 50(masked ladder but
                            // somethimes it has few hits), set bco_peak as -1
      {
        bco_peak = -1;
      }
      m_bcopeaks_fee[felix][fee] = bco_peak;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeBcoMapCdb()
{
  if (m_bcomap_cdb_file.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  CDBTTree *cdbttree = new CDBTTree(m_bcomap_cdb_file);

  int size = 0;
  std::vector<int> bco_temp_container;
  bco_temp_container.clear();
  for (int felix = 0; felix < 8; felix++)
  {
    if (m_FELIX_TARGET != -1 && m_FELIX_TARGET != felix)
    {
      continue;
    }
    std::string name_bco_peak = "h_InttCalib_BCOOffSet_INTT" + std::to_string(felix);
    m_bco_peak[felix] = new TH1I(name_bco_peak.c_str(), name_bco_peak.c_str(), 14, 0, 14);
    for (int fee = 0; fee < 14; fee++)
    {
      int bco = m_bcopeaks_fee[felix][fee];
      cdbttree->SetIntValue(size, "felix_server", felix);
      cdbttree->SetIntValue(size, "felix_channel", fee);
      cdbttree->SetIntValue(size, "bco_diff", bco);
      m_bco_peak[felix]->SetBinContent(fee + 1, bco);
      bco_temp_container.push_back(bco);
      ++size;
    }
    QAHistManagerDef::getHistoManager()->registerHisto(m_bco_peak[felix]);
  }

  cdbttree->SetSingleIntValue("size", size);
  std::pair<double, double> stats =
      CalculateStandardDeviation(bco_temp_container);
  m_bco_mean = stats.first;
  m_bco_stdDev = stats.second;
  int runmode = m_streaming ? 1 : 0;  //
  cdbttree->SetSingleDoubleValue("StdDev", m_bco_stdDev);
  cdbttree->SetSingleIntValue("runmode", runmode);
  cdbttree->SetSingleIntValue("events", m_evts);
  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();
  bco_temp_container.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::MakeBcoMapPng()
{
  // Canvas
  gStyle->SetOptStat(0);
  TCanvas *bco_cnvs = new TCanvas(  //
      "bco_cnvs", "bco_cnvs",       //
      1280, 720                     //
  );
  bco_cnvs->Draw();

  TH1 *bco_hist[112] = {};
  for (int i = 0; i < 16; ++i)
  {
    std::string tpadtitle = std::format("bco_pad_{:02d}", i);
    bco_cnvs->cd();
    TPad *bco_pdf_pad = new TPad(  //
        tpadtitle.c_str(), tpadtitle.c_str(),
        // NOLINTNEXTLINE(bugprone-integer-division)
        ((i % 4 + 0.0) / 4.0 * 0.9) + 0.0, ((3.0 - i / 4) / 4.0 * 0.9) + 0.1,
        // NOLINTNEXTLINE(bugprone-integer-division)
        ((i % 4 + 1.0) / 4.0 * 0.9) + 0.0, ((4.0 - i / 4) / 4.0 * 0.9) + 0.1);
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

    int felix = (i % 4) + (4 * (i / 8));
    int fee_index_start = (i / 4) % 2 ? 0 : 7;
    int fee_index_end = (i / 4) % 2 ? 7 : 14;
    double max = 0;
    for (int fee = fee_index_start; fee < fee_index_end; fee++)
    {
      int h = (felix * 14) + fee;
      bco_pdf_pad->cd();
      std::string htitle = std::format("h_InttCalib_bco_hist_{:03d}", h);
      // float min_bin = (m_streaming) ? -127.5 : -0.5;
      float min_bin = -0.5;
      float max_bin = 127.5;
      // int nbins = (m_streaming) ? 256 : 128;
      int nbins = 128;
      bco_hist[h] = new TH1D(              //
          htitle.c_str(), htitle.c_str(),  //
          nbins, min_bin, max_bin          //
                                           // 128, -0.5, 127.5                //
      );
      bco_hist[h]->SetTitle(std::format(";BCO Difference;intt{:01d} ({:01d} - {:02d})", felix, fee, fee).c_str());
      bco_hist[h]->GetYaxis()->CenterTitle();
      bco_hist[h]->GetYaxis()->SetTitleSize(0.12);
      bco_hist[h]->GetYaxis()->SetTitleOffset(0.6);
      bco_hist[h]->GetYaxis()->SetLabelSize(0.07);
      bco_hist[h]->GetXaxis()->SetTitleSize(0.10);
      bco_hist[h]->GetXaxis()->SetTitleOffset(0.6);
      bco_hist[h]->GetXaxis()->SetLabelSize(0.07);
      bco_hist[h]->SetLineColor(GetFeeColor(fee));
      bco_hist[h]->Draw("same");

      for (int bco = 0; bco < 128; ++bco)
      {
        bco_hist[h]->SetBinContent(bco + 1, m_bcorates_fee[felix][fee][bco]);
        max = std::max<double>(max, m_bcorates_fee[felix][fee][bco]);
      }
    }
    //	for (InttNameSpace::RawData_s raw = raw_begin; raw <= raw_end; ++raw)
    //	{
    //	  int h = (raw.felix_server) * 14 + raw.felix_channel;
    //	  bco_hist[h]->GetYaxis()->SetRangeUser(0.1, 10 * max);
    //	}
  }

  // Legend
  bco_cnvs->cd();
  TPad *legend_pad = new TPad("legend_pad", "legend_pad",  //
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
    legend_text.DrawText(
        0.6, (2.0 * fee + 1.0) / 14.0,
        std::format("FCh {:01d}, {:02d}", fee, (fee + 7)).c_str());

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
  TPad *caption_pad = new TPad("caption_pad", "caption_pad",  //
                               0.0, 0.0, 1.0, 0.1             //
  );
  caption_pad->SetFillStyle(4000);  // transparent
  caption_pad->Range(0.0, 0.0, 1.0, 1.0);
  caption_pad->Draw();

  caption_pad->cd();
  TText caption_text;
  caption_text.SetTextAlign(22);
  caption_text.SetTextSize(0.40);
  caption_text.SetTextColor(kBlack);
  caption_text.DrawText(
      0.5, 0.75,
      std::format("Run: {:08d} Events: {} BCO StdDev: {:.2f} BCO Offset: {:.2f}",
                  m_run_num, m_evts, m_bco_stdDev, m_bco_mean)
          .c_str());

  TText caption_tag;
  caption_tag.SetTextAlign(22);
  caption_tag.SetTextSize(0.40);
  if (m_evts < 1000)
  {
    caption_tag.SetTextColor(kRed);
    caption_tag.DrawText(0.5, 0.3, "LOW STATSTICS events < 1000!!");
  }
  else if (m_bco_stdDev != 0)
  {
    caption_tag.SetTextColor(kRed);
    caption_tag.DrawText(0.5, 0.3, "FEE MISALIGNED!!");
  }
  else
  {
    caption_tag.SetTextColor(kBlue);
    caption_tag.DrawText(0.5, 0.3, "GOOD");
  }

  std::map<int, std::vector<int>> maskedLadders;  // To store fee -> fc mappings

  // Iterate through the bco_hist array
  for (int i = 0; i < 112; ++i)
  {
    if (bco_hist[i]->GetBinContent(bco_hist[i]->GetMaximumBin()) < 5)
    {
      int sev = i / 14;
      int fee = i % 14;
      // Add fc to the corresponding fee in the map
      maskedLadders[sev].push_back(fee);
    }
  }

  // Now, prepare the text to add to the canvas
  std::ostringstream maskedText;
  maskedText << "Masked Ladders: ";

  // Iterate through the map and create the output text
  for (const auto &entry : maskedLadders)
  {
    maskedText << "INTT " << entry.first << "(FC ";  // Fee numbering starts at 1
    for (auto j = 0U; j < entry.second.size(); ++j)
    {
      maskedText << entry.second[j];
      if (j != entry.second.size() - 1)
      {
        maskedText << ",";  // Add commas between FC values
      }
    }
    maskedText << ") ";
  }

  // Print or draw the text on the canvas
  TText caption_masked;
  caption_masked.SetTextAlign(32);
  caption_masked.SetTextSize(0.20);
  caption_masked.DrawText(0.9, 0.3, maskedText.str().c_str());
  bco_cnvs->Update();
  bco_cnvs->Show();
  if (!m_bcomap_png_file.empty())
  {
    bco_cnvs->SaveAs(m_bcomap_png_file.c_str());
  }

  delete bco_cnvs;

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::SaveHitrates()
{
  TFile *file = TFile::Open(m_hotmap_png_file.c_str(), "RECREATE");
  if (!file)
  {
    std::cerr << "\n"
              << PHWHERE << "\n\tfile\n"
              << std::endl;
    return 1;
  }
  file->cd();
  TTree *tree = new TTree("hitrate_tree", "hitrate_tree");
  tree->SetDirectory(file);

  int pid{0};
  int fee{0};
  int chp{0};
  int chn{0};
  tree->Branch("pid", &pid);
  tree->Branch("fee", &fee);
  tree->Branch("chp", &chp);
  tree->Branch("chn", &chn);

  double hitrate;
  tree->Branch("hitrate", &hitrate);

  for (auto const &raw : InttNameSpace::AllRawDataChannels())
  {
    pid = raw.felix_server + 3001;
    fee = raw.felix_channel;
    chp = raw.chip;
    chn = raw.channel;

    hitrate = m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] / m_evts;
    tree->Fill();
  }

  tree->Write();
  file->Write();
  file->Close();

  return 0;
}

int InttCalib::LoadHitrates()
{
  TFile *file = TFile::Open("/sphenix/user/jbertaux/hitrates.root", "READ");
  if (!file)
  {
    std::cerr << "\n"
              << PHWHERE << "\n\tfile\n"
              << std::endl;
    return 1;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("hitrate_tree"));
  if (!tree)
  {
    std::cerr << "\n"
              << PHWHERE << "\n\ttree\n"
              << std::endl;
    return 1;
  }

  int pid{0};
  int fee{0};
  int chp{0};
  int chn{0};
  tree->SetBranchAddress("pid", &pid);
  tree->SetBranchAddress("fee", &fee);
  tree->SetBranchAddress("chp", &chp);
  tree->SetBranchAddress("chn", &chn);

  double hitrate;
  tree->SetBranchAddress("hitrate", &hitrate);

  for (Int_t n = 0, N = tree->GetEntriesFast(); n < N; ++n)
  {
    tree->GetEntry(n);
    InttNameSpace::RawData_s raw{
        .felix_server = pid + 3001,
        .felix_channel = fee,
        .chip = chp,
        .channel = chn,
    };
    m_hitmap[raw.felix_server][raw.felix_channel][raw.chip][raw.channel][128] = hitrate;
  }

  m_evts = 1.0;

  return 0;
}
int InttCalib::ConfigureHist_v3(TH1 *&hist, TF1 *&fit, double maxbin,
                                std::map<double, int> const &hitrate_map,
                                std::string const &name,
                                std::string const &title)
{
  size_t map_size = hitrate_map.size();
  double global_middle = maxbin / 5;  // Global middle value is maxbin / 5
  size_t mid_index = map_size / 2;
  double middle_key = 0.;

  auto it = hitrate_map.begin();
  std::advance(it, mid_index);
  middle_key = it->first;

  // Make hist
  delete hist;
  hist = new TH1D(                                        //
      (name + " hitrates").c_str(),                       //
      title.c_str(),                                      //
      100, std::next(hitrate_map.begin())->first, maxbin  //
  );

  TH1 *fit_hist = new TH1D(                               //
      (name + " fit_hitrates").c_str(),                   //
      (title + " (fit)").c_str(),                         //
      100, std::next(hitrate_map.begin())->first, maxbin  //
  );

  for (auto const &[hitrate, count] : hitrate_map)
  {
    for (int i = 0; i < count; ++i)
    {
      hist->Fill(hitrate);
      if (hitrate / global_middle > 0.7)
      {
        fit_hist->Fill(hitrate);
      }
    }
  }

  delete fit;
  fit = new TF1(                //
      (name + "_fit").c_str(),  //
      "gaus",                   //
      middle_key / 10, maxbin   //
  );

  if (Verbosity())
  {
    fit_hist->Fit(fit, "R");  // range, no-draw, log likelihood
  }
  else
  {
    fit_hist->Fit(fit, "QR");  // range, no-draw, log likelihood, quiet
  }

  std::cout << "Global middle: " << global_middle << std::endl;
  delete fit_hist;
  return 0;
}

int InttCalib::adjust_hitrate(InttNameSpace::Offline_s const &ofl,
                              double &hitrate)
{
  hitrate /= (ofl.ladder_z % 2) ? 2.0 : 1.6;  // normalize by sensor length
  if (ofl.layer < 5)
  {
    hitrate /= (10.005 / 7.4994);  // Inner = 7.4994, Outer = 10.005
  }

  // InttSurveyMap::val_t transform;
  // if (m_survey.GetStripTransform(ofl, transform))
  // {
  //   std::cout << PHWHERE << "\n"
  //             << "InttSurveyMap::GetStripTransform failed for\n"
  //             << ofl
  //             << std::endl;
  //   return 1;
  // }
  // Eigen::Vector3d v = transform.translation() - m_vertex;
  // Eigen::Vector3d n{transform.linear()(0, 1), transform.linear()(1, 1),
  // transform.linear()(2, 1)}; hitrate *= v.squaredNorm(); hitrate /= (-1.0 *
  // n.dot(v) / v.norm());
  return 0;
}

int InttCalib::GetIndex(InttNameSpace::RawData_s const &raw,
                        InttNameSpace::Offline_s const & /*ofl*/)
{
  // int index = 0;
  // index += (raw.felix_server) * 4;
  // index += (ofl.layer < 5) ? 0 : 2; // +2 for outer
  // index += (ofl.ladder_z % 2) ? 0 : 1; // +1 for type B
  // return index;
  return raw.felix_server;
}

int InttCalib::GetFeeIndex(InttNameSpace::RawData_s const &raw,
                           InttNameSpace::Offline_s const & /*ofl*/)
{
  int index = 0;
  index = (raw.felix_server) * 14 + raw.felix_channel;
  return index;
}
std::pair<double, double>
InttCalib::CalculateStandardDeviation(const std::vector<int> &data)
{
  if (data.empty())
  {
    return std::make_pair(-1, -1);
  }

  double sum = 0.0;
  int n_masked_ladder = 0;
  for (int i : data)
  {
    if (i == -1)
    {
      n_masked_ladder++;
      continue;
    }
    sum += i;
  }
  double mean = sum / static_cast<double>(data.size() - n_masked_ladder);
  double sumSquaredDiffs = 0.0;
  //  int count = 0;
  for (int i : data)
  {
    if (i == -1)
    {
      continue;  // do not include maksed ladder for std calculation
    }
    sumSquaredDiffs += (i - mean) * (i - mean);
    //    count++;
  }
  double stddev = std::sqrt(sumSquaredDiffs /
                            static_cast<double>(data.size() - n_masked_ladder));
  // return stddev;
  return std::make_pair(mean, stddev);
}

Color_t InttCalib::GetFeeColor(int fee)
{
  switch (fee % 7)
  {
  case 1:
    return kRed;
  case 2:
    return kGreen;
  case 3:
    return kBlue;
  case 4:
    return kOrange;
  case 5:
    return kMagenta;
  case 6:
    return kCyan;
  default:  // 0
    break;
  }
  return kBlack;
}
