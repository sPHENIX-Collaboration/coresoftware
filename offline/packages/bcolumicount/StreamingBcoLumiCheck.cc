#include "StreamingBcoLumiCheck.h"

//#include "BcoInfo.h"
#include "StreamingBcoInfo.h"
#include "StreamingLumiInfo.h"
//#include "BcoStreamingLumiInfov1.h"


#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllServer.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE


#include <iostream>

StreamingBcoLumiCheck::StreamingBcoLumiCheck(const std::string &name)
  : SubsysReco(name)
{
  return;
}

int StreamingBcoLumiCheck::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);

  return iret;
}

int StreamingBcoLumiCheck::InitRun(PHCompositeNode *topNode)
{
  StreamingLumiInfo *streaming_lumi_info = findNode::getClass<StreamingLumiInfo>(topNode, "STREAMINGLUMIINFO");
  if (streaming_lumi_info)
  {
    std::cout << " raw lumi : " << streaming_lumi_info->get_lumi_raw() << std::endl;
    std::cout << " live lumi : " << streaming_lumi_info->get_lumi_live() << std::endl;
    std::cout << " scaled lumi : " << streaming_lumi_info->get_lumi_scaled() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingBcoLumiCheck::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST Node is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingBcoLumiCheck::process_event(PHCompositeNode *topNode)
{
  StreamingBcoInfo *streaming_bco_info = findNode::getClass<StreamingBcoInfo>(topNode, "STREAMINGBCOINFO");
  if (streaming_bco_info)
  {
    if (Verbosity() > 1)
    {
      std::cout << "bco : " << streaming_bco_info->get_bco() << std::endl;
      std::cout << "usable bco tag : " << streaming_bco_info->get_usable_bco_tag() << std::endl;
      std::cout << "bco streaming window : (" <<  streaming_bco_info->get_bco_streaming_window().first << ", " << streaming_bco_info->get_bco_streaming_window().second << ")" << std::endl;
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}