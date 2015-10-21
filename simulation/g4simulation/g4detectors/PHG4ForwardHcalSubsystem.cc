#include "PHG4ForwardHcalSubsystem.h"
#include "PHG4ForwardHcalDetector.h"
#include "PHG4ForwardHcalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;


//_______________________________________________________________________
PHG4ForwardHcalSubsystem::PHG4ForwardHcalSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  active(1),
  detector_type(name),
  mappingfile_("")
{

}


//_______________________________________________________________________
int PHG4ForwardHcalSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4ForwardHcalDetector(topNode, Name());
  detector_->SetActive(active);
  detector_->OverlapCheck(overlapcheck);
  detector_->Verbosity(verbosity);
  detector_->SetTowerMappingFile( mappingfile_ );

  if (active)
    {
      // create hit output node
      ostringstream nodename;
      nodename <<  "G4HIT_" << detector_type;

      PHG4HitContainer* scintillator_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
      if (!scintillator_hits)
        {
          scintillator_hits = new PHG4HitContainer();
          PHIODataNode<PHObject> *hitNode = new PHIODataNode<PHObject>(scintillator_hits, nodename.str().c_str(), "PHObject");
          dstNode->addNode(hitNode);
        }

      ostringstream absnodename;
      absnodename << "G4HIT_ABSORBER_" << detector_type;

      PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
      if (!absorber_hits)
        {
          absorber_hits = new PHG4HitContainer();
          PHIODataNode<PHObject> *abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
          dstNode->addNode(abshitNode);
        }

      // create stepping action
      steppingAction_ = new PHG4ForwardHcalSteppingAction(detector_);
    }

  return 0;
}


//_______________________________________________________________________
int
PHG4ForwardHcalSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    {
      steppingAction_->SetInterfacePointers( topNode );
    }
  return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4ForwardHcalSubsystem::GetDetector( void ) const
{
  return detector_;
}


//_______________________________________________________________________
PHG4SteppingAction* PHG4ForwardHcalSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}
