#include "PHG4AccSubsystem.h"
#include "PHG4AccDetector.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4AccSteppingAction.h"
#include "PHG4CylinderEventAction.h"
#include <g4main/PHG4Utils.h>

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4AccSubsystem::PHG4AccSubsystem( const std::string &na, const int lyr):
  detector_( NULL ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  radius(100),
  length(100),
  xpos(0),
  ypos(0),
  zpos(0),
  lengthViaRapidityCoverage(true),
  TrackerThickness(100),
  material("Silicon"),
  _sciTilt(0),
  _sciWidth(0.1),
  _tungstenWidth(0.2),
  _unmag(1),
  _sciNum(-1),
  active(0),
  absorberactive(0),
  layer(lyr),
  detector_type(na),
  superdetector("NONE"),
  usedrawing(0)
{
  // put the layer into the name so we get unique names
  // for multiple SVX layers
  ostringstream nam;
  nam << na << "_" << lyr;
  Name(nam.str().c_str());
}

//_______________________________________________________________________
int PHG4AccSubsystem::InitRun( PHCompositeNode* topNode )
{
  // create hit list only for active layers
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));
  // create detector
  detector_ = new PHG4AccDetector(topNode, Name(), layer);
  detector_->SetRadius(radius);
  G4double detlength = length;
  if (lengthViaRapidityCoverage)
    {
      detlength =  PHG4Utils::GetLengthForRapidityCoverage(radius+TrackerThickness)*2;
    }
  detector_->SetLength(detlength);
  detector_->SetPosition(xpos, ypos, zpos);
  detector_->SetThickness(TrackerThickness);
  detector_->SetTilt(_sciTilt);
  detector_->SetScintWidth(_sciWidth);
  detector_->SetTungstenWidth(_tungstenWidth);
  detector_->SetUndulMag(_unmag);
  detector_->SetNumScint(_sciNum);
  detector_->SetMaterial(material);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  detector_->UseDrawing(usedrawing);
  if (active)
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
      PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str() );
      if ( !cylinder_hits )
        {
          dstNode->addNode( new PHIODataNode<PHObject>( cylinder_hits = new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));
        }
      cylinder_hits->AddLayer(layer);
      PHG4CylinderEventAction *evtac = new PHG4CylinderEventAction(topNode, nodename.str());
      if (absorberactive)
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
	  PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str() );
	  if ( !cylinder_hits )
	    {
	      dstNode->addNode( new PHIODataNode<PHObject>( cylinder_hits = new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));
	    }
	  cylinder_hits->AddLayer(layer);
          evtac->AddNode(nodename.str());
	}
      eventAction_ = evtac;
      steppingAction_ = new PHG4AccSteppingAction(detector_);
      steppingAction_->set_zmin(zpos-detlength/2.);
      steppingAction_->set_zmax(zpos + detlength/2.);
    }


  return 0;

}

//_______________________________________________________________________
int PHG4AccSubsystem::process_event( PHCompositeNode* topNode )
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
PHG4Detector* PHG4AccSubsystem::GetDetector( void ) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4AccSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}


void
PHG4AccSubsystem::Print(const std::string &what) const
{
  detector_->Print(what);
  return;
}
