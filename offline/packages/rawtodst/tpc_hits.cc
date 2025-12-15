
#include "tpc_hits.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitSetv1.h>
#include <trackbase/TrkrHitv2.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TH2.h>

//____________________________________________________________________________..
tpc_hits::tpc_hits(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "tpc_hits::tpc_hits(const std::string &name)" << std::endl;
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
tpc_hits::~tpc_hits()
{
  delete hm;
  std::cout << "tpc_hits::~tpc_hits() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int tpc_hits::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  hm = new Fun4AllHistoManager("HITHIST");

  _h_hit_XY = new TH2F("_h_hit_XY", "_h_hit_XY;X, [m];Y, [m]", 400, -800, 800, 400, -800, 800);

  hm->registerHisto(_h_hit_XY);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::InitRun(PHCompositeNode *topNode)
{
  // get dst node
  PHNodeIterator iter(topNode);
  auto *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "tpc_hits::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto *hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << "tpc_hits::InitRun - creating TRKR_HITSET." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto *trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    hitsetcontainer = new TrkrHitSetContainerv1;
    auto *newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }
  topNode->print();

  // we reset the BCO for the new run
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;

  // m_hits = new TrkrHitSetContainerv1();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::process_event(PHCompositeNode *topNode)
{
  // std::cout << "tpc_hits::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  //  load relevant nodes
  //  Get the TrkrHitSetContainer node
  auto *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  assert(_event);

  if (_event == nullptr)
  {
    std::cout << "tpc_hits::Process_Event - Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  // check all possible TPC packets that we need to analyze
  for (int ep = 0; ep < 2; ep++)
  {
    for (int sector = 0; sector <= 23; sector++)
    {
      const int packet = 4000 + sector * 10 + ep;

      // Reading packet
      Packet *p = _event->getPacket(packet);

      // Figure out which side
      int side = 0;
      if (sector > 11)
      {
        side = 1;
      }

      if (p)
      {
        std::cout << "tpc_hits:: Event getPacket: " << packet << "| Sector" << sector << "| EndPoint " << ep << std::endl;
      }
      else
      {
        continue;
      }

      int nr_of_waveforms = p->iValue(0, "NR_WF");

      for (auto &l : m_hitset)
      {
        l = new TrkrHitSetv1();

        int wf;
        for (wf = 0; wf < nr_of_waveforms; wf++)
        {
          int current_BCO = p->iValue(wf, "BCO") + rollover_value;

          if (starting_BCO < 0)
          {
            starting_BCO = current_BCO;
          }

          if (current_BCO < starting_BCO)  // we have a rollover
          {
            rollover_value += 0x100000;
            current_BCO = p->iValue(wf, "BCO") + rollover_value;
            starting_BCO = current_BCO;
            current_BCOBIN++;
          }
          int sampa_nr = p->iValue(wf, "SAMPAADDRESS");
          int channel = p->iValue(wf, "CHANNEL");

          int fee = p->iValue(wf, "FEE");
          int samples = p->iValue(wf, "SAMPLES");
          // clockwise FEE mapping
          // int FEE_map[26]={5, 6, 1, 3, 2, 12, 10, 11, 9, 8, 7, 1, 2, 4, 8, 7, 6, 5, 4, 3, 1, 3, 2, 4, 6, 5};
          int FEE_R[26] = {2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
          // conter clockwise FEE mapping (From Takao)
          int FEE_map[26] = {3, 2, 5, 3, 4, 0, 2, 1, 3, 4, 5, 7, 6, 2, 0, 1, 0, 1, 4, 5, 11, 9, 10, 8, 6, 7};
          int pads_per_sector[3] = {96, 128, 192};

          // setting the mapp of the FEE
          int feeM = FEE_map[fee] - 1;
          if (FEE_R[fee] == 2)
          {
            feeM += 6;
          }
          if (FEE_R[fee] == 3)
          {
            feeM += 14;
          }
          int layer = M.getLayer(feeM, channel);

          double R = M.getR(feeM, channel);
          double phi = M.getPhi(feeM, channel) + sector * M_PI / 6;

          TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(layer, sector, side);
          TrkrHitSetContainer::Iterator hitsetit = trkrhitsetcontainer->findOrAddHitSet(tpcHitSetKey);

          if (Verbosity())
          {
            int sampaAddress = p->iValue(wf, "SAMPAADDRESS");
            int sampaChannel = p->iValue(wf, "SAMPACHANNEL");
            int checksum = p->iValue(wf, "CHECKSUM");
            int checksumError = p->iValue(wf, "CHECKSUMERROR");
            std::cout << "tpc_hits::Process_Event Samples " << samples
                      << " Chn:" << channel
                      << " layer: " << layer
                      << " sampa: " << sampa_nr
                      << " fee: " << fee
                      << " Mapped fee: " << feeM
                      << " sampaAddress: " << sampaAddress
                      << " sampaChannel: " << sampaChannel
                      << " checksum: " << checksum
                      << " checksumError: " << checksumError
                      << " hitsetkey " << tpcHitSetKey
                      << " R = " << R
                      << " phi = " << phi
                      << std::endl;
          }
          for (int s = 0; s < samples; s++)
          {
            int pad = M.getPad(feeM, channel);
            int t = s + 2 * (current_BCO - starting_BCO);
            int adc = p->iValue(wf, s);
            // generate hit key
            TrkrDefs::hitkey hitkey = TpcDefs::genHitKey((unsigned int) pad + sector * pads_per_sector[FEE_R[sector] - 1], (unsigned int) t);
            // find existing hit, or create
            auto *hit = hitsetit->second->getHit(hitkey);

            // create hit, assign adc and insert in hitset
            if (!hit)
            {
              // create a new one
              hit = new TrkrHitv2();
              hit->setAdc(adc);
              hitsetit->second->addHitSpecificKey(hitkey, hit);
            }
            // else{
            //   hit->setAdc(adc);
            //   hitsetit->second->addHitSpecificKey(hitkey, hit);
            // }

            _h_hit_XY->Fill(R * cos(phi), R * sin(phi), adc);
          }
        }
      }

    }  // End of run over packets
  }    // End of ep loop
  // we skip the mapping to real pads at first. We just say
  // that we get 16 rows (segment R2) with 128 pads
  // so each FEE fills 2 rows. Not right, but one step at a time.

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int tpc_hits::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "tpc_hits::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  hm->dumpHistos(_filename, "RECREATE");

  return Fun4AllReturnCodes::EVENT_OK;
}
