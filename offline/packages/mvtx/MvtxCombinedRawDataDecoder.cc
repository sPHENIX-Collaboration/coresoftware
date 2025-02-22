/*!
 * \file MvtxCombinedRawDataDecoder.cc
 * \author Cameron Dean <cameron.dean@cern.ch>
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include "MvtxCombinedRawDataDecoder.h"

#include <fun4allraw/MvtxRawDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/MvtxEventInfov3.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContMvtxHelperv1.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <fun4all/Fun4AllServer.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/MvtxFeeIdInfov1.h>
#include <ffarawobjects/MvtxRawEvtHeaderv1.h>
#include <ffarawobjects/MvtxRawEvtHeaderv2.h>
#include <ffarawobjects/MvtxRawHitContainerv1.h>
#include <ffarawobjects/MvtxRawHitv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>
#include <phool/sphenix_constants.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>  // for accessing the MVTX hot pixel file from the CDB

#include <cassert>
#include <iterator>

namespace
{
  std::string MvtxHitSetHelperName("TRKR_MVTXHITSETHELPER");
}

//_________________________________________________________
MvtxCombinedRawDataDecoder::MvtxCombinedRawDataDecoder(const std::string &name)
  : SubsysReco(name)
{
}

//___________________________________________________________________________
void MvtxCombinedRawDataDecoder::CreateNodes(PHCompositeNode *topNode)
{
  // get dst node
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode =
      dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MvtxCombinedRawDataDecoder::InitRun - DST Node missing, "
                 "doing nothing."
              << std::endl;
    exit(1);
  }

  // create mvtx hitset helper container if needed
  mvtx_hit_set_helper =
      findNode::getClass<TrkrHitSetContMvtxHelper>(topNode, MvtxHitSetHelperName);
  if (!mvtx_hit_set_helper)
  {
    // find or create TRKR node
    auto *trkrNode = dynamic_cast<PHCompositeNode *>(
        iter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    // create container and add to the tree
    mvtx_hit_set_helper = new TrkrHitSetContMvtxHelperv1;
    auto *newNode = new PHIODataNode<PHObject>(mvtx_hit_set_helper, MvtxHitSetHelperName,
                                               "PHObject");
    trkrNode->addNode(newNode);
  }

  // create hitset container if needed
  hit_set_container =
      findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hit_set_container)
  {
    // find or create TRKR node
    auto *trkrNode = dynamic_cast<PHCompositeNode *>(
        iter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    // create container and add to the tree
    hit_set_container = new TrkrHitSetContainerv1;
    auto *newNode = new PHIODataNode<PHObject>(hit_set_container, "TRKR_HITSET",
                                               "PHObject");
    trkrNode->addNode(newNode);
  }

  // Check if MVTX event header already exists
  auto *mvtxNode = dynamic_cast<PHCompositeNode *>(
      iter.findFirst("PHCompositeNode", "MVTX"));
  if (!mvtxNode)
  {
    mvtxNode = new PHCompositeNode("MVTX");
    dstNode->addNode(mvtxNode);
  }

  mvtx_event_header =
      findNode::getClass<MvtxEventInfo>(mvtxNode, "MVTXEVENTHEADER");
  if (!mvtx_event_header)
  {
    mvtx_event_header = new MvtxEventInfov3();
    auto *newHeader = new PHIODataNode<PHObject>(
        mvtx_event_header, "MVTXEVENTHEADER", "PHObject");
    mvtxNode->addNode(newHeader);
  }

  mvtx_raw_event_header =
      findNode::getClass<MvtxRawEvtHeader>(topNode, m_MvtxRawEvtHeaderNodeName);

  mvtx_raw_hit_container =
      findNode::getClass<MvtxRawHitContainer>(topNode, m_MvtxRawHitNodeName);
}

//_____________________________________________________________________
void MvtxCombinedRawDataDecoder::GetNodes(PHCompositeNode *topNode)
{
  mvtx_raw_event_header =
      findNode::getClass<MvtxRawEvtHeader>(topNode, m_MvtxRawEvtHeaderNodeName);
  if (Verbosity() >= 3)
  {
    mvtx_raw_event_header->identify();
  }

  mvtx_raw_hit_container =
      findNode::getClass<MvtxRawHitContainer>(topNode, m_MvtxRawHitNodeName);
  if (Verbosity() >= 3)
  {
    mvtx_raw_hit_container->identify();
  }

  hit_set_container =
      findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  mvtx_hit_set_helper =
      findNode::getClass<TrkrHitSetContMvtxHelper>(topNode, MvtxHitSetHelperName);

  mvtx_event_header =
      findNode::getClass<MvtxEventInfo>(topNode, "MVTXEVENTHEADER");

  // Could we just get the first strobe BCO instead of setting this to 0?
  // Possible problem, what if the first BCO isn't the mean, then we'll shift tracker hit sets? Probably not a bad thing but depends on hit stripping
  //  uint64_t gl1rawhitbco = gl1 ? gl1->get_bco() : 0;
  auto *gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    gl1rawhitbco = gl1->lValue(0, "BCO");
  }
  else
  {
    auto *oldgl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
    if (oldgl1)
    {
      gl1rawhitbco = oldgl1->get_bco();
    }
  }
  if (gl1rawhitbco == 0 && (Verbosity() >= 4))
  {
    std::cout << PHWHERE << "Could not get gl1 raw hit" << std::endl;
  }
}

//_____________________________________________________________________
int MvtxCombinedRawDataDecoder::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MvtxCombinedRawDataDecoder::InitRun(PHCompositeNode *topNode)
{
  CreateNodes(topNode);

  if (!mvtx_raw_event_header || !mvtx_raw_hit_container)
  {
    Fun4AllServer::instance()->unregisterSubsystem(this);
    std::cout << PHWHERE << "::" << __func__ << ": Could not get \""
              << m_MvtxRawHitNodeName << " or " << m_MvtxRawEvtHeaderNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Unregistering subsystem and continuing on" << std::endl;

  }
  if (dynamic_cast<MvtxRawEvtHeaderv1 *>(mvtx_raw_event_header))
  {
    std::cout << PHWHERE <<  "MvtxCombinedRawDataDecoder::GetNodes() !!!WARNING!!! using obsolete MvtxRawEvtHeaderv1.";
    std::cout << " Unregistering MvtxCombinedRawDataDecoder SubsysReco." << std::endl;
    Fun4AllServer::instance()->unregisterSubsystem(this);
  }

  auto *rc = recoConsts::instance();
  int runNumber = rc->get_IntFlag("RUNNUMBER");

  if (m_readStrWidthFromDB)
  {
    m_strobeWidth = MvtxRawDefs::getStrobeLength(runNumber);
  }
  if (std::isnan(m_strobeWidth))
  {
    std::cout << "MvtxCombinedRawDataDecoder::InitRun - strobe width is undefined for this run, defaulting to 89 mus" << std::endl;
    m_strobeWidth = 89;
  }
  if (m_strobeWidth < 1)
  {
    runMvtxTriggered(true);
  }
  //! save MVTX strobe width
  rc->set_FloatFlag("MvtxStrobeWidth", m_strobeWidth);

  // Load the hot pixel map from the CDB
  if (m_doOfflineMasking)
  {
    m_hot_pixel_mask = new MvtxPixelMask();
    m_hot_pixel_mask->load_from_CDB();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int MvtxCombinedRawDataDecoder::process_event(PHCompositeNode *topNode)
{
  GetNodes(topNode);

  // get the last 40 bits by bit shifting left then right to match
  // to the mvtx bco
  auto lbshift = gl1rawhitbco << 24U;
  auto gl1bco = lbshift >> 24U;

  // std::vector<std::pair<uint64_t, uint32_t> > strobe_bc_pairs;
  // std::set<uint64_t> l1BCOs = mvtx_raw_event_header->getMvtxLvL1BCO();
  // auto mvtxbco = *l1BCOs.begin();
  // if (gl1rawhitbco && (Verbosity() > 2))
  // {
  //   std::cout << "MVTX header BCO " << mvtxbco << " and GL1 BCO " << gl1bco
  //             << std::endl;
  // }

  for (const auto &L1 : mvtx_raw_event_header->getMvtxLvL1BCO())
  {
    mvtx_event_header->add_L1_BCO(L1);
  }
  auto nMvtxFeeIdInfo = mvtx_raw_event_header->get_nFeeIdInfo();
  for (size_t i{}; i < nMvtxFeeIdInfo; ++i)
  {
    auto *feeIdInfo = mvtx_raw_event_header->get_feeIdInfo(i);
    mvtx_event_header->add_strobe_BCO(feeIdInfo->get_bco());
  }
  if (Verbosity() >= 3)
  {
    mvtx_raw_hit_container->identify();
  }
  const auto &strobe_list = mvtx_event_header->get_strobe_BCOs();
  auto it_strb_bco_zero = strobe_list.upper_bound(gl1bco);
  auto str_wGL1_idx = std::distance(strobe_list.cbegin(), it_strb_bco_zero) - 1;

  uint64_t hit_strobe = -1;  // Initialise to -1 for debugging
  uint8_t layer = 0;
  uint8_t stave = 0;
  uint8_t chip = 0;
  uint16_t row = 0;
  uint16_t col = 0;

  for (unsigned int i = 0; i < mvtx_raw_hit_container->get_nhits(); i++)
  {
    mvtx_rawhit = mvtx_raw_hit_container->get_hit(i);
    hit_strobe = mvtx_rawhit->get_bco();
    layer = mvtx_rawhit->get_layer_id();
    stave = mvtx_rawhit->get_stave_id();
    chip = mvtx_rawhit->get_chip_id();
    row = mvtx_rawhit->get_row();
    col = mvtx_rawhit->get_col();

    int strobe_index = -20;
    const auto it = strobe_list.find(hit_strobe);
    if (it == strobe_list.cend())
    {
      std::cout << "Warning: hit strobe BCO " << hit_strobe << " is not found in evet combined strobe list:" << std::endl;
      for (const auto &strobe : strobe_list)
      {
        std::cout << "0x" << std::hex << strobe << std::dec << std::endl;
      }
    }
    else
    {
      strobe_index = static_cast<int>(std::distance(strobe_list.cbegin(), it) - str_wGL1_idx);
    }

    // int bcodiff = gl1rawhitbco ? strobe - gl1bco : 0;
    //   double timeElapsed = bcodiff * 0.1065;  // 106 ns rhic clock
    //   int index = m_mvtx_is_triggered ? 0 : std::ceil(timeElapsed / m_strobeWidth);

    if (strobe_index < -16 || strobe_index > 15)
    {
      std::cout << "Strobe index: " << strobe_index << " out of range" << std::endl;
      continue;  // Index is out of the 5-bit signed range
    }

    if (Verbosity() >= 10)
    {
      mvtx_rawhit->identify();
    }

    const TrkrDefs::hitsetkey hitsetkey =
        MvtxDefs::genHitSetKey(layer, stave, chip, strobe_index);
    if (!hitsetkey)
    {
      continue;
    }

    mvtx_hit_set_helper->addHitSetKey(strobe_index, hitsetkey);

    // get matching hitset
    const auto hitset_it = hit_set_container->findOrAddHitSet(hitsetkey);

    // generate hit key
    const TrkrDefs::hitkey hitkey = MvtxDefs::genHitKey(col, row);

    // find existing hit, or create
    auto *hit = hitset_it->second->getHit(hitkey);
    if (hit)
    {
      if (Verbosity() > 1)
      {
        std::cout << PHWHERE << "::" << __func__
                  << " - duplicated hit, hitsetkey: " << hitsetkey
                  << " hitkey: " << hitkey << std::endl;
      }
      continue;
    }

    if (m_doOfflineMasking)
    {
      if (!m_hot_pixel_mask->is_masked(mvtx_rawhit))
      {  // Check if the pixel is masked
        hit = new TrkrHitv2;
        hitset_it->second->addHitSpecificKey(hitkey, hit);
      }
    }
    else
    {
      hit = new TrkrHitv2;
      hitset_it->second->addHitSpecificKey(hitkey, hit);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MvtxCombinedRawDataDecoder::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


// void MvtxCombinedRawDataDecoder::removeDuplicates(
//     std::vector<std::pair<uint64_t, uint32_t> > &v)
// {
//   auto end = v.end();
//   for (auto it = v.begin(); it != end; ++it)
//   {
//     end = remove(it + 1, end, *it);
//   }
//   v.erase(end, v.end());
// }
