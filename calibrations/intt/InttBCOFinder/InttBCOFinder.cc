#include "InttBCOFinder.h"

#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

InttBCOFinder::InttBCOFinder(const std::string &name, const std::string &filename, const std::string &filename2, int nevents)
  : SubsysReco(name)
  , nevents_(nevents)
  , outfname_(filename)
  , cdbname_(filename2)
{
}

// Destructor
InttBCOFinder::~InttBCOFinder()
{
  for (int i = 0; i < 8; i++)
  {
    delete h2_bco_ladder_[i];
    delete h2_bco_ladder_cut_[i];
  }
}
/**
 * Initialize the module and prepare looping over events
 */
int InttBCOFinder::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in InttBCOFinder" << std::endl;
  }
  for (int j = 0; j < 8; j++)
  {
    h2_bco_ladder_[j] = new TH2D((boost::format("h2_bco_felix_%d") % j).str().c_str(), (boost::format("h2_bco_felix_%d") % j).str().c_str(), 14, 0, 14, 128, 0, 128);
    h2_bco_ladder_cut_[j] = new TH2D((boost::format("h2_bco_felix_cut%d") % j).str().c_str(), (boost::format("h2_bco_felix_cut%d") % j).str().c_str(), 14, 0, 14, 128, 0, 128);
  }
  return 0;
}

int InttBCOFinder::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning InitRun in InttBCOFinder" << std::endl;
  }

  if (!topNode)
  {
    std::cout << "InttBCOFinder::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;

    return 1;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttBCOFinder::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "InttBCOFinder::Beginning process_event in InttBCOFinder.cc" << std::endl;
  }

  if (ievent_ >= nevents_)
  {
    if (Verbosity() > 5)
    {
      std::cout << "InttBCOFinder::Last event is processed." << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
    // return Fun4AllReturnCodes::ABORTRUN;
  }

  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "InttBCOFinder::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit *intthit = inttcont->get_hit(i);
    // if (!intthit)
    //   continue;
    uint64_t bco_full = intthit->get_bco();
    int bco = intthit->get_FPHX_BCO();
    int felixnumber = intthit->get_packetid() - 3001;
    int felixchannel = intthit->get_fee();
    // std::cout << bco << " " << felixnumber << std::endl;
    int bco_diff = (bco_full & 0x7FU) - bco;

    if (bco_diff < 0)
    {
      h2_bco_ladder_[felixnumber]->Fill(felixchannel, bco_diff + 128);
    }
    else
    {
      h2_bco_ladder_[felixnumber]->Fill(felixchannel, bco_diff);
    }
  }
  ievent_++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttBCOFinder::End(PHCompositeNode * /*topNode*/)
{
  FindBCOPeak();
  if (WriteCDBTTree_)
  {
    cdbttree_ = new CDBTTree(cdbname_);
    int size = 0;
    if (Verbosity() > 1)
    {
      std::cout << "InttBCOFinder::Creating CDBTTree..." << std::endl;
    }
    for (int felix_server = 0; felix_server < 8; felix_server++)
    {
      for (int felix_channel = 0; felix_channel < 14; felix_channel++)
      {
        for (int bco_diff = 0; bco_diff < 128; bco_diff++)
        {
          if (h2_bco_ladder_cut_[felix_server]->GetBinContent(felix_channel + 1, bco_diff + 1) > 0)
          {
            cdbttree_->SetIntValue(size, "felix_server", felix_server);
            cdbttree_->SetIntValue(size, "felix_channel", felix_channel);
            cdbttree_->SetIntValue(size, "bco_diff", bco_diff);
            size++;
          }
        }
      }
    }

    cdbttree_->SetSingleIntValue("size", size);
    cdbttree_->Commit();
    cdbttree_->CommitSingle();
    cdbttree_->WriteCDBTTree();
  }

  delete cdbttree_;

  if (Verbosity() > 1)
  {
    std::cout << "InttBCOFinder::Processing InttBCOFinder done" << std::endl;
  }
  if (WriteQAFile_)
  {
    if (Verbosity() > 1)
    {
      std::cout << "InttBCOFinder::Writing histograms of BCO distribution" << std::endl;
      std::cout << "InttBCOFinder::File path : " << outfname_ << std::endl;
    }
    outFile_ = TFile::Open(outfname_.c_str(), "RECREATE");
    if (outFile_ != nullptr)
    {
      outFile_->cd();
      for (int j = 0; j < 8; j++)
      {
        h2_bco_ladder_[j]->Write();
        h2_bco_ladder_cut_[j]->Write();
      }
      outFile_->Write();
      outFile_->Close();
    }
  }
  return 0;
}

void InttBCOFinder::FindBCOPeak()
{
  for (int felixnumber = 0; felixnumber < 8; felixnumber++)
  {
    int numXBins = h2_bco_ladder_[felixnumber]->GetXaxis()->GetNbins();
    int numYBins = h2_bco_ladder_[felixnumber]->GetYaxis()->GetNbins();

    // Find the maximum bin half ladder by half ladder//
    for (int binX = 1; binX <= numXBins; binX++)
    {
      double maxXValue = -1.0;
      int maxYBin = -1;
      for (int binY = 1; binY <= numYBins; binY++)
      {
        double binContent = h2_bco_ladder_[felixnumber]->GetBinContent(binX, binY);
        if (binContent > maxXValue)
        {
          maxXValue = binContent;
          maxYBin = binY;
        }
      }
      h2_bco_ladder_cut_[felixnumber]->SetBinContent(binX, maxYBin, maxXValue);  // Fill the peak position in the 2D histogram( it will be used for BCO cut)
      // Check the closest bin content (peak-1, peak+1)
      int BinY1 = maxYBin - 1;
      if (maxYBin == 1)
      {
        BinY1 = 128;
      }
      int BinY2 = maxYBin + 1;
      if (maxYBin == 128)
      {
        BinY2 = 1;
      }
      double ClosestBinContent1 = h2_bco_ladder_[felixnumber]->GetBinContent(binX, BinY1);
      double ClosestBinContent2 = h2_bco_ladder_[felixnumber]->GetBinContent(binX, BinY2);
      h2_bco_ladder_cut_[felixnumber]->SetBinContent(binX, BinY1, ClosestBinContent1);
      h2_bco_ladder_cut_[felixnumber]->SetBinContent(binX, BinY2, ClosestBinContent2);
    }
    h2_bco_ladder_[felixnumber]->SetTitle((boost::format("felix %d evt %d") % felixnumber % ievent_).str().c_str());
    h2_bco_ladder_[felixnumber]->SetXTitle("module");
    h2_bco_ladder_[felixnumber]->SetYTitle("BCO_Full - BCO");
    h2_bco_ladder_cut_[felixnumber]->SetTitle((boost::format("felix %d evt %d cut") % felixnumber % ievent_).str().c_str());
    h2_bco_ladder_cut_[felixnumber]->SetXTitle("module");
    h2_bco_ladder_cut_[felixnumber]->SetYTitle("BCO_Full - BCO");
  }
}
