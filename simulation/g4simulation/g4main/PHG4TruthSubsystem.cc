#include "PHG4TruthSubsystem.h"

#include "PHG4Particle.h"                // for PHG4Particle
#include "PHG4TruthEventAction.h"
#include "PHG4TruthTrackingAction.h"
#include "PHG4TruthInfoContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHObject.h>              // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                 // for PHWHERE

#include <cassert>
#include <cstdlib>                      // for exit
#include <iostream>
#include <set>                           // for _Rb_tree_iterator, set, _Rb_...
#include <utility>                       // for pair

class PHG4EventAction;
class PHG4TrackingAction;

using namespace std;

//_______________________________________________________________________
PHG4TruthSubsystem::PHG4TruthSubsystem(const string& name)
  : PHG4Subsystem(name)
  , m_EventAction(nullptr)
  , m_TrackingAction(nullptr)
  , m_SaveOnlyEmbededFlag(false)
{
}

//_______________________________________________________________________
int PHG4TruthSubsystem::InitRun(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create truth information container
  PHG4TruthInfoContainer* truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthInfoList)
  {
    truthInfoList = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(truthInfoList, "G4TruthInfo", "PHObject"));
  }

  // event action
  m_EventAction = new PHG4TruthEventAction();

  // create tracking action
  m_TrackingAction = new PHG4TruthTrackingAction(m_EventAction);

  return 0;
}

//_______________________________________________________________________
int PHG4TruthSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to event action so that it gets
  // relevant nodes needed internally
  if (m_EventAction)
  {
    m_EventAction->SetInterfacePointers(topNode);
  }
  else
  {
    cout << PHWHERE << " No EventAction registered" << endl;
    exit(1);
  }

  if (m_TrackingAction)
  {
    m_TrackingAction->SetInterfacePointers(topNode);
  }
  else
  {
    cout << PHWHERE << " No TrackingAction registered" << endl;
    exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TruthSubsystem::process_after_geant(PHCompositeNode* topNode)
{
  if (m_SaveOnlyEmbededFlag)
  {
    if (Verbosity() > 1)
    {
      cout << __PRETTY_FUNCTION__ << " - INFO - only save the G4 truth information that is associated with the embedded particle" << endl;
    }

    PHG4TruthInfoContainer* truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    assert(truthInfoList);

    set<int> savevtxlist;

    // remove particle that is not embedd associated
    PHG4TruthInfoContainer::Range truth_range = truthInfoList->GetParticleRange();
    PHG4TruthInfoContainer::Iterator truthiter = truth_range.first;
    while (truthiter != truth_range.second)
    {
      const int primary_id = (truthiter->second)->get_primary_id();
      if (truthInfoList->isEmbeded(primary_id) <= 0)
      {
        // not a embed associated particle

        truthInfoList->delete_particle(truthiter++);
      }
      else
      {
        // save vertex id for primary particle which leaves no hit
        // in active area
        savevtxlist.insert((truthiter->second)->get_vtx_id());
        ++truthiter;
      }
    }

    // remove vertex that is not embedd associated
    PHG4TruthInfoContainer::VtxRange vtxrange = truthInfoList->GetVtxRange();
    PHG4TruthInfoContainer::VtxIterator vtxiter = vtxrange.first;
    while (vtxiter != vtxrange.second)
    {
      if (savevtxlist.find(vtxiter->first) == savevtxlist.end())
      {
        truthInfoList->delete_vtx(vtxiter++);
      }
      else
      {
        ++vtxiter;
      }
    }
  }

  return 0;
}

int PHG4TruthSubsystem::ResetEvent(PHCompositeNode* topNode)
{
  m_TrackingAction->ResetEvent(topNode);
  m_EventAction->ResetEvent(topNode);
  return 0;
}

//_______________________________________________________________________
PHG4EventAction* PHG4TruthSubsystem::GetEventAction(void) const
{
  return m_EventAction;
}

//_______________________________________________________________________

PHG4TrackingAction*
PHG4TruthSubsystem::GetTrackingAction(void) const
{
  return m_TrackingAction;
}
