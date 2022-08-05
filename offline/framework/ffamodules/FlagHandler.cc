#include "FlagHandler.h"

#include <ffaobjects/FlagSave.h>
#include <ffaobjects/FlagSavev1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/recoConsts.h>

//____________________________________________________________________________..
FlagHandler::FlagHandler(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int FlagHandler::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  FlagSave *flgsv = findNode::getClass<FlagSave>(runNode, "Flags");
  if (!flgsv)
  {
    flgsv = new FlagSavev1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(flgsv, "Flags", "PHObject");
    runNode->addNode(newNode);
  }
  else
  {
    recoConsts *rc = recoConsts::instance();
    flgsv->PutFlagsBack(rc, false);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int FlagHandler::End(PHCompositeNode *topNode)
{
  FlagSave *flagsave = findNode::getClass<FlagSave>(topNode, "Flags");
  if (flagsave)
  {
    recoConsts *rc = recoConsts::instance();
    flagsave->FillFromPHFlag(rc, true);
    flagsave->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void FlagHandler::Print(const std::string & /* what */) const
{
  recoConsts *rc = recoConsts::instance();
  rc->Print();
}
