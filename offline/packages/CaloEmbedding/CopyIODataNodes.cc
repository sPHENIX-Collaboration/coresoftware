#include "CopyIODataNodes.h"

#include <globalvertex/GlobalVertexMap.h>

#include <calotrigger/MinimumBiasInfo.h>

#include <centrality/CentralityInfo.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>
#include <ffaobjects/SyncObject.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

//____________________________________________________________________________..
CopyIODataNodes::CopyIODataNodes(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int CopyIODataNodes::InitRun(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  if (m_CopyRunHeaderFlag)
  {
    CopyRunHeader(topNode, se->topNode());
  }
  if (m_CopyEventHeaderFlag)
  {
    CreateEventHeader(topNode, se->topNode());
  }
  if (m_CopyCentralityInfoFlag)
  {
    CreateCentralityInfo(topNode, se->topNode());
  }
  if (m_CopyGlobalVertexMapFlag)
  {
    CreateGlobalVertexMap(topNode, se->topNode());
  }
  if (m_CopyMinimumBiasInfoFlag)
  {
    CreateMinimumBiasInfo(topNode, se->topNode());
  }
  if (m_CopySyncObjectFlag)
  {
    CreateSyncObject(topNode, se->topNode());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CopyIODataNodes::process_event(PHCompositeNode *topNode)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  if (m_CopyEventHeaderFlag)
  {
    CopyEventHeader(topNode, se->topNode());
  }
  if (m_CopyCentralityInfoFlag)
  {
    CopyCentralityInfo(topNode, se->topNode());
  }
  if (m_CopyGlobalVertexMapFlag)
  {
    CopyGlobalVertexMap(topNode, se->topNode());
  }
  if (m_CopyMinimumBiasInfoFlag)
  {
    CopyMinimumBiasInfo(topNode, se->topNode());
  }
  if (m_CopySyncObjectFlag)
  {
    CopySyncObject(topNode, se->topNode());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CopyIODataNodes::CopyRunHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  RunHeader *from_runheader = findNode::getClass<RunHeader>(from_topNode, "RunHeader");
  if (!from_runheader)
  {
    std::cout << "Could not locate RunHeader on " << from_topNode->getName() << std::endl;
    m_CopyRunHeaderFlag = false;
    return;
  }
  RunHeader *to_runheader = findNode::getClass<RunHeader>(to_topNode, "RunHeader");
  if (!to_runheader)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    if (!runNode)
    {
      runNode = new PHCompositeNode("RUN");
      to_topNode->addNode(runNode);
    }
    to_runheader = dynamic_cast<RunHeader *>(from_runheader->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_runheader, "RunHeader", "PHObject");
    runNode->addNode(newNode);
    if (Verbosity() > 0)
    {
      std::cout << "From RunHeader identify()" << std::endl;
      from_runheader->identify();
      std::cout << "To RunHeader identify()" << std::endl;
      to_runheader->identify();
    }
  }
  return;
}

void CopyIODataNodes::CreateCentralityInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  CentralityInfo *from_centralityinfo = findNode::getClass<CentralityInfo>(from_topNode, "CentralityInfo");
  if (!from_centralityinfo)
  {
    std::cout << "Could not locate CentralityInfo on " << from_topNode->getName() << std::endl;
    m_CopyCentralityInfoFlag = false;
    return;
  }
  CentralityInfo *to_centralityinfo = findNode::getClass<CentralityInfo>(to_topNode, "CentralityInfo");
  if (!to_centralityinfo)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      to_topNode->addNode(dstNode);
    }
    to_centralityinfo = dynamic_cast<CentralityInfo *>(from_centralityinfo->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_centralityinfo, "CentralityInfo", "PHObject");
    dstNode->addNode(newNode);
  }
  return;
}

void CopyIODataNodes::CopyCentralityInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  CentralityInfo *from_centralityinfo = findNode::getClass<CentralityInfo>(from_topNode, "CentralityInfo");
  CentralityInfo *to_centralityinfo = findNode::getClass<CentralityInfo>(to_topNode, "CentralityInfo");
  from_centralityinfo->CopyTo(to_centralityinfo);
  if (Verbosity() > 0)
  {
    std::cout << "From CentralityInfo identify()" << std::endl;
    from_centralityinfo->identify();
    std::cout << "To CentralityInfo identify()" << std::endl;
    to_centralityinfo->identify();
  }

  return;
}

void CopyIODataNodes::CreateEventHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  EventHeader *from_eventheader = findNode::getClass<EventHeader>(from_topNode, "EventHeader");
  if (!from_eventheader)
  {
    std::cout << "Could not locate EventHeader on " << from_topNode->getName() << std::endl;
    m_CopyEventHeaderFlag = false;
    return;
  }
  EventHeader *to_eventheader = findNode::getClass<EventHeader>(to_topNode, "EventHeader");
  if (!to_eventheader)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      to_topNode->addNode(dstNode);
    }
    to_eventheader = dynamic_cast<EventHeader *>(from_eventheader->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_eventheader, "EventHeader", "PHObject");
    dstNode->addNode(newNode);
  }
  return;
}

void CopyIODataNodes::CopyEventHeader(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  EventHeader *from_eventheader = findNode::getClass<EventHeader>(from_topNode, "EventHeader");
  EventHeader *to_eventheader = findNode::getClass<EventHeader>(to_topNode, "EventHeader");
  from_eventheader->CopyTo(to_eventheader);
  if (Verbosity() > 0)
  {
    std::cout << "From EventHeader identify()" << std::endl;
    from_eventheader->identify();
    std::cout << "To EventHeader identify()" << std::endl;
    to_eventheader->identify();
  }

  return;
}

void CopyIODataNodes::CreateGlobalVertexMap(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  GlobalVertexMap *from_globalvertexmap = findNode::getClass<GlobalVertexMap>(from_topNode, "GlobalVertexMap");
  if (!from_globalvertexmap)
  {
    std::cout << "Could not locate GlobalVertexMap on " << from_topNode->getName() << std::endl;
    m_CopyGlobalVertexMapFlag = false;
    return;
  }
  GlobalVertexMap *to_globalvertexmap = findNode::getClass<GlobalVertexMap>(to_topNode, "GlobalVertexMap");
  if (!to_globalvertexmap)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      to_topNode->addNode(dstNode);
    }
    to_globalvertexmap = dynamic_cast<GlobalVertexMap *>(from_globalvertexmap->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_globalvertexmap, "GlobalVertexMap", "PHObject");
    dstNode->addNode(newNode);
  }
  return;
}

void CopyIODataNodes::CopyGlobalVertexMap(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  GlobalVertexMap *from_globalvertexmap = findNode::getClass<GlobalVertexMap>(from_topNode, "GlobalVertexMap");
  GlobalVertexMap *to_globalvertexmap = findNode::getClass<GlobalVertexMap>(to_topNode, "GlobalVertexMap");
  from_globalvertexmap->CopyTo(to_globalvertexmap);
  if (Verbosity() > 0)
  {
    std::cout << "From GlobalVertexMap identify()" << std::endl;
    from_globalvertexmap->identify();
    std::cout << "To GlobalVertexMap identify()" << std::endl;
    to_globalvertexmap->identify();
  }

  return;
}

void CopyIODataNodes::CreateMinimumBiasInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  MinimumBiasInfo *from_minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(from_topNode, "MinimumBiasInfo");
  if (!from_minimumbiasinfo)
  {
    std::cout << "Could not locate MinimumBiasInfo on " << from_topNode->getName() << std::endl;
    m_CopyMinimumBiasInfoFlag = false;
    return;
  }
  MinimumBiasInfo *to_minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(to_topNode, "MinimumBiasInfo");
  if (!to_minimumbiasinfo)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      to_topNode->addNode(dstNode);
    }
    to_minimumbiasinfo = dynamic_cast<MinimumBiasInfo *>(from_minimumbiasinfo->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_minimumbiasinfo, "MinimumBiasInfo", "PHObject");
    dstNode->addNode(newNode);
  }
  return;
}

void CopyIODataNodes::CopyMinimumBiasInfo(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  MinimumBiasInfo *from_minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(from_topNode, "MinimumBiasInfo");
  MinimumBiasInfo *to_minimumbiasinfo = findNode::getClass<MinimumBiasInfo>(to_topNode, "MinimumBiasInfo");
  from_minimumbiasinfo->CopyTo(to_minimumbiasinfo);
  if (Verbosity() > 0)
  {
    std::cout << "From MinimumBiasInfo identify()" << std::endl;
    from_minimumbiasinfo->identify();
    std::cout << "To MinimumBiasInfo identify()" << std::endl;
    to_minimumbiasinfo->identify();
  }

  return;
}

void CopyIODataNodes::CreateSyncObject(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  SyncObject *from_syncobject = findNode::getClass<SyncObject>(from_topNode, "Sync");
  if (!from_syncobject)
  {
    std::cout << "Could not locate SyncObject on " << from_topNode->getName() << std::endl;
    m_CopySyncObjectFlag = false;
    return;
  }
  SyncObject *to_syncobject = findNode::getClass<SyncObject>(to_topNode, "Sync");
  if (!to_syncobject)
  {
    PHNodeIterator iter(to_topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      to_topNode->addNode(dstNode);
    }
    to_syncobject = dynamic_cast<SyncObject *>(from_syncobject->CloneMe());
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(to_syncobject, "SyncObject", "PHObject");
    dstNode->addNode(newNode);
  }
  return;
}

void CopyIODataNodes::CopySyncObject(PHCompositeNode *from_topNode, PHCompositeNode *to_topNode)
{
  SyncObject *from_syncobject = findNode::getClass<SyncObject>(from_topNode, "Sync");
  SyncObject *to_syncobject = findNode::getClass<SyncObject>(to_topNode, "Sync");
  to_syncobject = from_syncobject;
  if (Verbosity() > 0)
  {
    std::cout << "From Syncobject identify()" << std::endl;
    from_syncobject->identify();
    std::cout << "To Syncobject identify()" << std::endl;
    to_syncobject->identify();
  }
  return;
}
