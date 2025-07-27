#include "EventFlowReco.h"

#include "EventFlowInfo.h"
#include "EventFlowInfov1.h"

#include "EventFlowInfoMap.h"
#include "EventFlowInfoMapv1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <epd/EpdGeom.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h> // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h> // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h> // for PHWHERE
#include <phool/recoConsts.h>


EventFlowReco::EventFlowReco(const std::string &name) 
  : SubsysReco(name) 
{ 
}

int EventFlowReco::InitRun(PHCompositeNode *topNode) 
{

  if(Verbosity() > 1) {
    std::cout << "EventFlowReco::InitRun -- entered" << std::endl;
  }
  return CreateNodes(topNode);
}

int EventFlowReco::process_event(PHCompositeNode *topNode) 
{

  if (Verbosity() > 1) {
    std::cout << "EventFlowReco::process_event -- entered" << std::endl;
  }
  if(!topNode) {
    std::cout << "EventFlowReco::process_event -- topNode is null, exiting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventFlowReco::CreateNodes(PHCompositeNode *topNode) {
  
  PHNodeIterator iter(topNode);

  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto globalNode = dynamic_cast<PHCompositeNode *>( iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode) {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  auto eps = findNode::getClass<EventFlowInfoMap>(topNode, "EventFlowInfoMap");
  if (!eps)
  {
    eps = new EventFlowInfoMapv1();
    auto EpMapNode = new PHIODataNode<PHObject>(eps, "EventFlowInfoMap", "PHObject");
    globalNode->addNode(EpMapNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}



int EventFlowReco::End(PHCompositeNode * /*topNode*/) {
  
  if (Verbosity() > 1) {
    std::cout << "EventFlowReco::End -- entered" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
