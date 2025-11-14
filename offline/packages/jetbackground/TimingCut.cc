#include "TimingCut.h"

#include <mbd/MbdOut.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <cmath>
#include <iostream>  // for basic_ostream, operator<<
#include <map>       // for _Rb_tree_iterator, opera...
#include <utility>   // for pair
#include <vector>    // for vector
//____________________________________________________________________________..
TimingCut::TimingCut(const std::string &jetNodeName, const std::string &name, const bool doAbort)
  : SubsysReco(name)
  , _doAbort(doAbort)
  , _jetNodeName(jetNodeName)
  , _cutParams(name)
{
  SetDefaultParams();
}

//____________________________________________________________________________..
int TimingCut::Init(PHCompositeNode *topNode)
{
  if(CreateNodeTree(topNode))
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TimingCut::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  if (!parNode)
  {
    std::cout << "No RUN node found; cannot create PHParameters for storing cut results!";
    return 1;
  }

  _cutParams.SaveToNodeTree(parNode, "TimingCutParams");
  return 0;
}

//____________________________________________________________________________..
int TimingCut::process_event(PHCompositeNode *topNode)
{
  JetContainer *jets = findNode::getClass<JetContainer>(topNode, _jetNodeName);
  if (!jets)
  {
    if (Verbosity() > 0 && !_missingInfoWarningPrinted)
    {
      std::cout << "Missing jets; abort event. Further warnings will be suppressed." << std::endl;
    }
    _missingInfoWarningPrinted = true;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  float maxJetpT = 0;
  float subJetpT = 0;
  float maxJett = std::numeric_limits<float>::quiet_NaN();
  float subJett = std::numeric_limits<float>::quiet_NaN();
  
  if (jets)
  {
    int tocheck = jets->size();
    if (Verbosity() > 2)
    {
      std::cout << "Found " << tocheck << " jets to check..." << std::endl;
    }
    for (int i = 0; i < tocheck; ++i)
    {
      float jetpT = 0;
      float jett = std::numeric_limits<float>::quiet_NaN();
      Jet *jet = jets->get_jet(i);
      if (jet)
      {
        jetpT = jet->get_pt();
	jett = jet->get_property(Jet::PROPERTY::prop_t);
      }
      else
      {
        continue;
      }
      if (jetpT > maxJetpT)
      {
        if (maxJetpT)
        {
          subJetpT = maxJetpT;
          subJett = maxJett;
        }
        maxJetpT = jetpT;
        maxJett = jett;
      }
      else if (jetpT > subJetpT)
      {
        subJetpT = jetpT;
        subJett = jett;
      }
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "No jet node!" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  bool passDeltat = Pass_Delta_t(maxJett, subJett);
  bool passLeadt = Pass_Lead_t(maxJett);

  MbdOut * mbdout = static_cast<MbdOut*>(findNode::getClass<MbdOut>(topNode,"MbdOut"));
  float m_mbd_t0 = std::numeric_limits<float>::quiet_NaN();
  float m_mbd_ts = std::numeric_limits<float>::quiet_NaN();
  float m_mbd_tn = std::numeric_limits<float>::quiet_NaN();
  if(mbdout)
    {
      m_mbd_t0 = mbdout->get_t0();
      m_mbd_ts = mbdout->get_time(0); // south side
      m_mbd_tn = mbdout->get_time(1); // north side
    }

  float mbd_time = std::numeric_limits<float>::quiet_NaN();
  if(!std::isnan(m_mbd_t0))
    {
      mbd_time = m_mbd_t0;
    }
  else if(!std::isnan(m_mbd_tn))
    {
      mbd_time = m_mbd_tn;
    }
  else if(!std::isnan(m_mbd_ts))
    {
      mbd_time = m_mbd_ts;
    }
  
  bool passMbdt = false;
  if(!std::isnan(mbd_time))
    {
      passMbdt = Pass_Mbd_dt(maxJett, mbd_time);
    }

  bool failAnyCut = !passDeltat || !passLeadt || !passMbdt;
    
  if (failAnyCut && _doAbort)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  _cutParams.set_int_param("passLeadtCut",passLeadt);
  _cutParams.set_int_param("passDeltatCut",passDeltat);
  _cutParams.set_int_param("passMbdDtCut",passMbdt);
  _cutParams.set_int_param("failAnyTimeCut", failAnyCut);
  _cutParams.set_double_param("maxJett",maxJett);
  _cutParams.set_double_param("subJett",subJett);
  _cutParams.set_double_param("mbd_time",mbd_time);
  _cutParams.UpdateNodeTree(parNode, "TimingCutParams");

  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int TimingCut::ResetEvent(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "TimingCut::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TimingCut::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "TimingCut::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TimingCut::Reset(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "TimingCut::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void TimingCut::Print(const std::string &what) const
{
  std::cout << "TimingCut::Print(const std::string &what) const Printing info for " << what << std::endl;
}
