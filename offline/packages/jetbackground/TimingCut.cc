#include "TimingCut.h"

#include <calobase/RawTowerDefs.h>  // for CalorimeterId, encode_to...
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/Vertex.h>  // for Vertex

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
TimingCut::TimingCut(const std::string &jetNodeName, const std::string &name, const int debug, const bool doAbort, GlobalVertex::VTXTYPE vtxtype)
  : SubsysReco(name)
  , _doAbort(doAbort)
  , _name(name)
  , _debug(debug)
  , _jetNodeName(jetNodeName)
  , _vtxtype(vtxtype)
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
    if (_debug > 0 && !_missingInfoWarningPrinted)
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
    if (_debug > 2)
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
    if (_debug > 0)
    {
      std::cout << "No jet node!" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  bool failsDeltat = Fails_Delta_t(maxJett, subJett);
  bool failsLeadt = Fails_Lead_t(maxJett);
  bool failsMbdt = Fails_Mbd_dt(maxJett, subJett);

  bool failsAnyCut = failsDeltat || failsLeadt || failsMbdt;
    
  if (failsAnyCut && _doAbort)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  _cutParams.set_int_param("failsLeadtCut",failsLeadt);
  _cutParams.set_int_param("failsDeltatCut",failsDeltat);
  _cutParams.set_int_param("failsMbdDtCut",failsMbdt);
  _cutParams.set_int_param("failsAnyTimeCut", failsAnyCut);
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
