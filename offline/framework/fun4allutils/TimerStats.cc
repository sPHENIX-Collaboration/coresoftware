#include "TimerStats.h"

#include <cdbobjects/CDBTTree.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/getClass.h>

#include <map>  // for _Rb_tree_iterator

TimerStats::TimerStats(const std::string &name)
  : SubsysReco(name)
{
}

int TimerStats::InitRun(PHCompositeNode * /*topNode*/)
{
  delete cdbttree; // make cppcheck happy, deleting a null ptr
  cdbttree = new CDBTTree(outfilename);
  return Fun4AllReturnCodes::EVENT_OK;
}

int TimerStats::process_event(PHCompositeNode *topNode)
{
  iev++;
  Fun4AllServer *se = Fun4AllServer::instance();
  EventHeader *evtheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (evtheader)
  {
    iev = evtheader->get_EvtSequence();
  }
  else
  {
    iev++;
  }
  for (auto iter = se->timer_begin(); iter != se->timer_end(); ++iter)
  {
    cdbttree->SetFloatValue(iev, iter->first, iter->second.elapsed());
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TimerStats::End(PHCompositeNode * /*topNode*/)
{
  cdbttree->Commit();
  cdbttree->WriteCDBTTree();
  delete cdbttree;
  return Fun4AllReturnCodes::EVENT_OK;
}
