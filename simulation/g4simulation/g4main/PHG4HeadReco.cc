#include "PHG4HeadReco.h"

#include "PHG4EventHeaderv1.h"

#include <ffaobjects/FlagSavev1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <phool/phool.h>

#include <TSystem.h>

#include <HepMC/GenEvent.h>

#include <iostream>

using namespace std;

PHG4HeadReco::PHG4HeadReco(const std::string &name):
  SubsysReco(name),
  evtseq(0)
{}	   

int
PHG4HeadReco::Init(PHCompositeNode *topNode)
{
  enum {DSTNODE, RUNNODE, LAST}; // leave LAST at end - it is used for loops
  // first test if neccessary nodes have been created, if not bail out
  const char *NName[] = {
    "DST",
    "RUN"};

  PHNodeIterator iter(topNode);
  PHCompositeNode *outNode[LAST];
  for (int i=0; i<LAST; i++)
    {
      outNode[i] = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",NName[i]));
      if (!outNode[i])
	{
	  cout << PHWHERE << NName[i] << " node is missing, no point in continuing exiting now" << endl;
	  gSystem->Exit(1);
  // just to make scan-build happy which does not know that gSystem->Exit() terminates
	  exit(1);
	}
    }	 
  PHG4EventHeader *eventheader = new PHG4EventHeaderv1();
  PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(eventheader,"EventHeader","PHObject");
  outNode[DSTNODE]->addNode(newNode);

  FlagSave *flgsv = new FlagSavev1();
  PHIODataNode<PHObject> *FlagSaveNode =
    new PHIODataNode<PHObject>(flgsv, "Flags", "PHObject");
  outNode[RUNNODE]->addNode(FlagSaveNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4HeadReco::process_event(PHCompositeNode *topNode)
{
  HepMC::GenEvent *hepmcevt = findNode::getClass<HepMC::GenEvent>(topNode, "HEPMC");
  PHG4EventHeader *evtheader = findNode::getClass<PHG4EventHeader>(topNode,"EventHeader");
  evtseq++;
  if (hepmcevt)
    {
      evtseq = hepmcevt->event_number();
      HepMC::HeavyIon *hi = hepmcevt->heavy_ion();
      if (hi)
	{
	  evtheader->set_ImpactParameter(hi->impact_parameter());
	  evtheader->set_EventPlaneAngle(hi->event_plane_angle());
	}
    }
  evtheader->set_EvtSequence(evtseq);
  if (Verbosity() > 0)
    {
      evtheader->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
