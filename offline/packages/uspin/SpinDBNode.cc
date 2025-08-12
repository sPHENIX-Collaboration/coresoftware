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
#include "SpinDBContentv1.h"
#include "SpinDBOutput.h"

SpinDBNode::SpinDBNode(const std::string &name)
  : SubsysReco(name)
{
}

int SpinDBNode::InitRun(PHCompositeNode *topNode)
{
  RunHeader *runHeader = findNode::getClass<RunHeader>(topNode, "RunHeader");

  if (!runHeader)
  {
    std::cout << PHWHERE << ":: RunHeader node missing! Skipping run XXX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter(topNode);

  PHCompositeNode *lowerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("RUN");
    topNode->addNode(lowerNode);
  }

  SpinDBContent *spin_cont = findNode::getClass<SpinDBContent>(topNode, "SpinDBContent");
  if (!spin_cont)
  {
    int runnumber = runHeader->get_RunNumber();
    spin_cont = new SpinDBContentv1;
    SpinDBOutput spin_out("phnxrc");
    spin_out.StoreDBContent(runnumber, runnumber);
    spin_out.GetDBContentStore(spin_cont, runnumber);
    PHIODataNode<PHObject> *spindbcontentNode = new PHIODataNode<PHObject>(spin_cont, "SpinDBContent", "PHObject");
    lowerNode->addNode(spindbcontentNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
