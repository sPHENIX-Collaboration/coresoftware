#include "PHG4InnerHcalSubsystem.h"
#include "PHG4InnerHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4InnerHcalSteppingAction.h"
#include "PHG4InnerHcalParameters.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4InnerHcalSubsystem::PHG4InnerHcalSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_(NULL),
  steppingAction_( NULL ),
  eventAction_(NULL),
  layer(lyr),
  detector_type(name),
  superdetector("NONE")
{

  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  params = new PHG4InnerHcalParameters();
}

//_______________________________________________________________________
int PHG4InnerHcalSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4InnerHcalDetector(topNode, params, Name());
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  set<string> nodes;
  if (params->IsActive())
    {
      ostringstream nodename;
      if (superdetector != "NONE")
	{
	  nodename <<  "G4HIT_" << superdetector;
	}
      else
	{
	  nodename <<  "G4HIT_" << detector_type << "_" << layer;
	}
      nodes.insert(nodename.str());
      if (params->IsAbsorberactive())
	{
	  nodename.str("");
	  if (superdetector != "NONE")
	    {
	      nodename <<  "G4HIT_ABSORBER_" << superdetector;
	    }
	  else
	    {
	      nodename <<  "G4HIT_ABSORBER_" << detector_type << "_" << layer;
	    }
          nodes.insert(nodename.str());
	}
      BOOST_FOREACH(string node, nodes)
	{
	  PHG4HitContainer* g4_hits =  findNode::getClass<PHG4HitContainer>( topNode , node.c_str());
	  if ( !g4_hits )
	    {
	      g4_hits = new PHG4HitContainer();
	      dstNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));

	    }
	  if (! eventAction_)
	    {
	      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, node);
	    }
	  else
	    {
	      PHG4EventActionClearZeroEdep *evtact = dynamic_cast<PHG4EventActionClearZeroEdep *>(eventAction_);

	      evtact->AddNode(node);
	    }
	}

      // create stepping action
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, params);

    }
  else
    {
      // if this is a black hole it does not have to be active
      if (params->IsBlackHole())
	{
	  steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, params);
	}
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4InnerHcalSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    {
      steppingAction_->SetInterfacePointers( topNode );
    }
  return 0;
}


void
PHG4InnerHcalSubsystem::Print(const string &what) const
{
  cout << "Inner Hcal Parameters: " << endl;
  params->print();
  if (detector_)
    {
      detector_->Print(what);
    }
  return;
}

//_______________________________________________________________________
PHG4Detector* PHG4InnerHcalSubsystem::GetDetector( void ) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4InnerHcalSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}

PHG4InnerHcalParameters *
PHG4InnerHcalSubsystem::GetParameters()
{
  return params;
}

void
PHG4InnerHcalSubsystem::SetActive(const int i)
{
  params->SetActive(i);
}

void
PHG4InnerHcalSubsystem::SetAbsorberActive(const int i)
{
  params->SetAbsorberactive(i);
}

void
PHG4InnerHcalSubsystem::BlackHole(const int i)
{
  params->BlackHole(i);
}

void
PHG4InnerHcalSubsystem::SetTiltViaNcross(const int ncross)
{
  if (ncross == 0)
    {
      cout << "Invalid number of crossings: " << ncross
	   << " how do you expect me to calculate a tilt angle for this????"
	   << endl
	   << "If you want a 0 degree tilt angle, just use SetTiltAngle(0)"
	   << endl
	   << "I refuse to continue this!" << endl;
      exit(1);
    }
  params->set_ncross(ncross);
}

