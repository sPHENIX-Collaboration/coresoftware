#include "SyncReco.h"

#include <ffaobjects/SyncObject.h>
#include <ffaobjects/SyncObjectv1.h>
#include <ffaobjects/SyncDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHIODataNode.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>


#include <iostream>

using namespace std;

static const char *SyncNodeName = "Sync";

SyncReco::SyncReco(const string &name): SubsysReco(name)
{
  return ;
}

int SyncReco::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  return iret;
}

int SyncReco::InitRun(PHCompositeNode *topNode)
{
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  recoConsts *rc = recoConsts::instance();
  syncobject->RunNumber(rc->get_IntFlag("RUNNUMBER"));
  return Fun4AllReturnCodes::EVENT_OK;
}

int SyncReco::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      cout << PHWHERE << " DST Node is missing doing nothing" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, SyncNodeName);
  if (!syncobject)
    {
      syncobject = new SyncObjectv1();
      PHIODataNode <PHObject> *SyncObjectNode = new PHIODataNode <PHObject>(syncobject, SyncNodeName, "PHObject"); // contains PHObject
      dstNode->addNode(SyncObjectNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SyncReco::process_event(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, SyncNodeName);
  if (!syncobject)
    {
      cout << PHWHERE << " No Synchronisation Object, no parallel reading of multiple inputs" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  syncobject->EventCounter(se->EventCounter());
  syncobject->EventNumber(se->EventNumber());
  syncobject->RunNumber(se->RunNumber());
  syncobject->SegmentNumber(se->SegmentNumber());

  if (Verbosity() > 0)
    {
      syncobject->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

