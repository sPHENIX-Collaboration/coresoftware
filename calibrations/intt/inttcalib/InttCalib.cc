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
#include <TH1D.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>

#include <boost/format.hpp>

#include <cmath>
#include <iostream>
#include <limits>

InttCalib::InttCalib(const std::string& name)
  : SubsysReco(name)
{
}

int InttCalib::InitRun(PHCompositeNode* /*unused*/)
{
  m_evts = 0;
  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    for (int bco = 0; bco < 129; ++bco)
    {
      m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] = 0;
    }
  }

  if (m_survey.LoadFromCDB("InttSurveyMap"))
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not load 'InttSurveyMap' from CDB\n"
              << "\tExiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  if (m_feemap.LoadFromCDB("InttFeeMap"))
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not load 'InttFeeMap' from CDB\n"
              << "\tExiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::process_event(PHCompositeNode* top_node)
{
  InttRawHitContainer* intt_raw_hit_container = findNode::getClass<InttRawHitContainer>(top_node, "INTTRAWHIT");
  if (!intt_raw_hit_container)
  {
    std::cout << PHWHERE << "\n"
              << "\tCould not get 'INTTRAWHIT' from node tree\n"
              << "\tExiting" << std::endl;
    gSystem->Exit(1);
    exit(1);
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (size_t n = 0, N = intt_raw_hit_container->get_nhits(); n < N; ++n)
  {
    InttRawHit* intt_raw_hit = intt_raw_hit_container->get_hit(n);
    if (!intt_raw_hit)
    {
      std::cout << PHWHERE << "\n"
                << "\tInttRawHit is nullptr but in range of InttRawHitContainer::get_nhits\n"
                << "\tExiting" << std::endl;
      gSystem->Exit(1);
      exit(1);
    }

    InttMap::RawData_s raw{
        //
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
  m_run_num = run_number;

  ConfigureHotMap();
  MakeHotMapCdb();
  MakeHotMapPng();

  ConfigureBcoMap();
  MakeBcoMapCdb();
  MakeBcoMapPng();

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCalib::ConfigureHotMap()
{
  m_hitrates.clear();
  for (InttMap::RawData_s raw = InttMap::RawDataBegin; raw != InttMap::RawDataEnd; ++raw)
  {
    double hitrate = (double) m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / (double) m_evts;
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
    return Fun4AllReturnCodes::ABORTEVENT;
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
    return Fun4AllReturnCodes::ABORTEVENT;
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
    double hitrate = (double) m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / (double) m_evts;
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
    double hitrate = (double) m_hitmap[raw.pid - 3001][raw.fee][raw.chp][raw.chn][128] / (double) m_evts;
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
    return Fun4AllReturnCodes::ABORTEVENT;
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
    return Fun4AllReturnCodes::ABORTEVENT;
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

void InttCalib::Debug() const
{
  InttSurveyMap::val_t transform;
  if (m_survey.GetStripTransform(InttMap::OfflineBegin, transform))
  {
    std::cout << PHWHERE << "\n"
              << "InttSurveyMap::GetStripTransform failed\n"
              << std::endl;
    return;
  }
  std::cout << "translation:" << std::endl;
  std::cout << transform.translation() << std::endl;

  std::cout << "vertex:" << std::endl;
  std::cout << m_vertex << std::endl;

  std::cout << "difference:" << std::endl;
  Eigen::Vector3d v = transform.translation() - m_vertex;
  std::cout << v << std::endl;

  std::cout << "normal:" << std::endl;
  Eigen::Vector3d n{transform.linear()(0, 1), transform.linear()(1, 1), transform.linear()(2, 1)};
  std::cout << n << std::endl;
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
