#include "InttHitMap.h"

#include <intt/InttMap.h>

#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <boost/format.hpp>

#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// Constructor
InttHitMap::InttHitMap(const std::string &name, const std::string &filename, int nevents)
  : SubsysReco(name)
  , nevents_(nevents)
  , outfname_(filename)
{
}

// Destructor
InttHitMap::~InttHitMap()
{
  for (auto &i : h2_AllMap_)
  {
    for (auto &j : i)
    {
      delete j;
    }
  }
}

/**
 * Initialize the module and prepare looping over events
 */
int InttHitMap::Init(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in InttHitMap" << std::endl;
  }
  outFile_ = new TFile(outfname_.c_str(), "RECREATE");
  for (int i = 0; i < 8; i++)
  {
    for (int j = 0; j < 14; j++)
    {
      h2_AllMap_[i][j] = new TH2D((boost::format("HitMap_%d_%d") % i % j).str().c_str(), (boost::format("Normalized_HitMap_%d_%d") % i % j).str().c_str(), 128, 0, 128, 27, 0, 27);
    }
  }
  return 0;
}

int InttHitMap::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning InitRun in InttHitMap" << std::endl;
  }

  if (!topNode)
  {
    std::cout << "InttHitMap::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;

    return -1;
  }

  if (nevents_ == 0)
  {
    std::cout << "InttHitMap::Maximun # of event is not specified." << std::endl;
    std::cout << "InttHitMap::Total # of events will be used." << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttHitMap::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "InttHitMap::Beginning process_event in InttHitMap" << std::endl;
  }
  //  nevents_ = 1000;
  if (ievent_ >= nevents_ && nevents_ != 0)
  {
    std::cout << "Last event is processed." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
    // return Fun4AllReturnCodes::ABORTRUN;
  }
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "InttHitMap::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }
  memset(IsCloneHit_, false, sizeof(IsCloneHit_));

  InttMap::RawData_s raw;
  InttMap::Online_s onl;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit *intthit = inttcont->get_hit(i);
    // if (!intthit)
    //   continue;
    raw.pid = intthit->get_packetid();
    raw.fee = intthit->get_fee();
    raw.chp = (intthit->get_chip_id() + 25) % 26;
    raw.chn = intthit->get_channel_id();

    int in_chip = raw.chp;
    int in_chan = raw.chn;
    int felixnumber = raw.pid - 3001;
    int in_module = raw.fee;

    if (fee_map.Convert(onl, raw))
    {
      std::cerr << "failed to convert raw:\n"
                << raw << "\n(continuing)" << std::endl;
      continue;
    }
    int layer = onl.lyr;
    int in_barrel = layer / 2;

    uint64_t bco_full = intthit->get_bco();
    int bco = intthit->get_FPHX_BCO();

    //////////////////// BCO Filtering  ////////////////////////////////
    if (isBCOcutON_ && (inBCOFile_ != nullptr))
    {
      if (!isBCOPeak(felixnumber, in_module, bco, bco_full))
      {
        continue;
      }
    }
    ///////////////////////////////////////////////////////////////////
    if (IsCloneHit_[felixnumber][in_module][in_chip][in_chan])
    {
      continue;
    }
    else
    {
      IsCloneHit_[felixnumber][in_module][in_chip][in_chan] = true;
    }
    FillHitMap(felixnumber, in_module, in_barrel, in_chip, in_chan);
  }
  ievent_++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttHitMap::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Processing InttHitMap done" << std::endl;
  }
  if (outFile_ != nullptr)
  {
    outFile_->Write();
    outFile_->Close();
  }
  return 0;
}

int InttHitMap::SetFeeMapFile(const std::string &feemapfile)
{
  if (fee_map.LoadFromFile(feemapfile))
  {
    std::cerr << "InttHitMap::failed to load fee map" << std::endl;
    return 1;
  }
  return 0;
}

int InttHitMap::SetBCOFile(const std::string &bcofile)
{
  if (!isBCOcutON_ && Verbosity() > 5)
  {
    std::cout << "InttHitMap::BCO cut option is OFF. isBCOcutON_ == false)" << std::endl;
    return 0;
  }
  inBCOFile_ = TFile::Open(bcofile.c_str());
  if (inBCOFile_ == nullptr)
  {
    std::cout << "InttHitMap::BCO file is not sucessfully loaded." << std::endl;
    return 0;
  }

  for (int i = 0; i < 8; i++)
  {
    inBCOFile_->GetObject((boost::format("h2_bco_felix_cut%d") % i).str().c_str(), h2_bco_cut_[i]);
  }
  return 1;
}

bool InttHitMap::isBCOPeak(int felix, int ladder, int bco, uint64_t bcofull)
{
  int bco_diff = (bcofull & 0x7FU) - bco;
  if (bco_diff < 0)
  {
    bco_diff = bco_diff + 128;
  }
  // Main part for BCO cut.
  // Load the BCO histogram and apply the BCO cut
  // Default : Hits belongs to [peak-1,peak+1] (3BCO) region will survice after BCO cut
  // To change the BCO region, h2_bco_cut_[felix] needs to be modified
  if (h2_bco_cut_[felix]->GetBinContent(ladder + 1, bco_diff + 1) != 0)
  {
    return true;
  }
  return false;
}

bool InttHitMap::FillHitMap(int in_felix, int in_module, int in_barrel, int in_chip, int in_chan)
{
  double norm_factor = 0.;
  if (in_chip < 5 || (in_chip > 12 && in_chip < 18))  // Condition for Type B
  {
    norm_factor = (nevents_);
    if (in_barrel == 1)  // Inner barrel(= 0) : 7.1888[cm], Outer barrel(= 1) : 9.680[cm]
    {
      norm_factor = norm_factor / (9.680 / 7.1888);
    }
    h2_AllMap_[in_felix][in_module]->Fill(in_chan, in_chip, 1. / norm_factor);
  }
  else
  {
    norm_factor = (nevents_) / 1.25;  // Normailzed by chip type
    if (in_barrel == 1)               // Inner barrel(= 0) : 7.1888[cm], Outer barrel(= 1) : 9.680[cm]
    {
      norm_factor = norm_factor / (9.680 / 7.1888);
    }
    h2_AllMap_[in_felix][in_module]->Fill(in_chan, in_chip, 1. / (norm_factor));
  }
  return true;
}
