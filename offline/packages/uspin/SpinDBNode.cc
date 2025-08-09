#include "SpinDBNode.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h> 
#include <phool/getClass.h>

#include <ffaobjects/RunHeader.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/SubsysReco.h>

#include "SpinDBContent.h"
#include "SpinDBOutput.h"

SpinDBNode::SpinDBNode(const std::string &name):
 SubsysReco(name)
{
}

SpinDBNode::~SpinDBNode()
{
}

int SpinDBNode::Init(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SpinDBNode::InitRun(PHCompositeNode *topNode)
{
  RunHeader* runHeader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  
  if (!runHeader)
  {
    std::cout << PHWHERE << ":: RunHeader node missing! Skipping run XXX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  int runnumber = runHeader->get_RunNumber();
  SpinDBContent spin_cont;
  SpinDBOutput spin_out("phnxrc");
  spin_out.StoreDBContent(runnumber, runnumber);
  spin_out.GetDBContentStore(spin_cont, runnumber);

  PHNodeIterator iter(topNode);

  PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("RUN");
    topNode->addNode(lowerNode);
  }

  PHIODataNode<SpinDBContent>* spindbcontentNode = new PHIODataNode<SpinDBContent>(&spin_cont, "SpinDBContent");
  lowerNode->addNode(spindbcontentNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SpinDBNode::process_event(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SpinDBNode::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

