#include "Fun4AllMemoryHistograms.h"

#include "Fun4AllHistoManager.h"
#include "Fun4AllUtils.h"

#include <TH1.h>
#include <TSystem.h>

#include <iostream>
#include <utility>  // for pair, make_pair

Fun4AllMemoryHistograms *Fun4AllMemoryHistograms::mInstance {nullptr};

Fun4AllMemoryHistograms::Fun4AllMemoryHistograms()
  : Fun4AllBase("Fun4AllMemoryHistograms")
{
}

void Fun4AllMemoryHistograms::FillHistos()
{
  if (! m_Enabled)
  {
    return;
  }
  if (! m_HistosCreated)
  {
    CreateHistos();
  }

  Fun4AllUtils::memoryvals mymem;
  Fun4AllUtils::GetMemory(mymem);
  if (m_EventCount >= m_NumChannels)
  {
    increase_TH1_channelcount();
  }
  int EventCountPlus1 = m_EventCount+1; // histos start at channel 1 ....
  m_histovector[0]->SetBinContent(EventCountPlus1,mymem.arena);
  m_histovector[1]->SetBinContent(EventCountPlus1,mymem.uordblks);
  m_histovector[2]->SetBinContent(EventCountPlus1,mymem.fordblks);
  m_histovector[3]->SetBinContent(EventCountPlus1,mymem.hblkhd);
  m_histovector[4]->SetBinContent(EventCountPlus1,mymem.VmHWM);
  m_histovector[5]->SetBinContent(EventCountPlus1,mymem.VmRSS);
  m_EventCount++;
}

void Fun4AllMemoryHistograms::SaveHistos(const std::string &fname)
{
  if (! m_Enabled)
  {
    return;
  }
  Fun4AllHistoManager *hm = new Fun4AllHistoManager("MEMHISTOS");
  for (auto iter : m_histovector)
  {
    hm->registerHisto(iter);
  }
  hm->dumpHistos(fname);
  return;
}

void Fun4AllMemoryHistograms::increase_TH1_channelcount()
{
  m_NumChannels += 10000;
  float m_Range = m_NumChannels + 0.5;
  for (TH1* histo : m_histovector)
  {
    const int old_nbins = histo->GetNbinsX();
    std::vector<double> old_content(old_nbins + 2);
// TH1::SetBins wipes the content, so we need to save it
    for (int bin = 0; bin <= old_nbins + 1; ++bin)
    {
      old_content[bin] = histo->GetBinContent(bin);
    }

    histo->SetBins(m_NumChannels, 0.5, m_Range);  // set number of bins which wipes the old histo content

// and then restore it
    for (int bin = 0; bin <= old_nbins + 1; ++bin)
    {
      histo->SetBinContent(bin, old_content[bin]);
    }
  }
  return;
}

void Fun4AllMemoryHistograms::CreateHistos()
{
  float m_Range = m_NumChannels + 0.5;
  TH1 *arena = new TH1D("arena","arena size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(arena);
  TH1 *uordblks = new TH1D("uordblks","uordblks size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(uordblks);
  TH1 *fordblks = new TH1D("fordblks","fordblks size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(fordblks);
  TH1 *hblkhd = new TH1D("hblkhd","hblkhd size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(hblkhd);
  TH1 *VmHWM = new TH1D("VmHWM","VmHWM size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(VmHWM);
  TH1 *VmRSS = new TH1D("VmRSS","VmRSS size",m_NumChannels,0.5,m_Range);
  m_histovector.push_back(VmRSS);
  m_HistosCreated = true;
}
