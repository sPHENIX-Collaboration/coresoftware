#include "PHG4CylinderSubsystem.h"
#include "PHG4CylinderDetector.h"
#include "PHG4CylinderGeomv1.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4FlushStepTrackingAction.h"
#include "PHG4CylinderRegionSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4CylinderSubsystem::PHG4CylinderSubsystem( const std::string &na, const int lyr):
  detector_( NULL ),
  steppingAction_( NULL ),
  trackingAction_(NULL),
  eventAction_(NULL),
  radius(100),
  length(100),
  xpos(0),
  ypos(0),
  zpos(0),
  lengthViaRapidityCoverage(true),
  TrackerThickness(100),
  material("G4_Galactic"), // default - almost nothing
  active(0),
  reduced(false),
  layer(lyr),
  blackhole(0),
  detector_type(na),
  superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple SVX layers
  ostringstream nam;
  nam << na << "_" << lyr;
  Name(nam.str().c_str());
}

//_______________________________________________________________________
int PHG4CylinderSubsystem::InitRun( PHCompositeNode* topNode )
{
  // create hit list only for active layers
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));
  // create detector
  detector_ = new PHG4CylinderDetector(topNode, Name(), layer);
  detector_->SetRadius(radius);
  G4double detlength = length;
  if (lengthViaRapidityCoverage)
    {
      detlength =  PHG4Utils::GetLengthForRapidityCoverage(radius+TrackerThickness)*2;
    }
  detector_->SetLength(detlength);
  detector_->SetPosition(xpos, ypos, zpos);
  detector_->SetThickness(TrackerThickness);
  detector_->SetMaterial(material);
  detector_->SetActive(active);
  detector_->BlackHole(blackhole);
  detector_->SetReducedTruthInfo(reduced);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  if (active)
    {
      ostringstream nodename;
      ostringstream geonode;
      if (superdetector != "NONE")
        {
          nodename <<  "G4HIT_" << superdetector;
          geonode << "CYLINDERGEOM_" << superdetector;
        }
      else
        {
          nodename <<  "G4HIT_" << detector_type << "_" << layer;
          geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
        }
      PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str() );
      if ( !cylinder_hits )
        {
          dstNode->addNode( new PHIODataNode<PHObject>( cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));
        }
      cylinder_hits->AddLayer(layer);
      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonode.str().c_str());
      if (!geo)
        {
          geo = new PHG4CylinderGeomContainer();
          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
          runNode->addNode(newNode);
        }
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv1(radius, zpos-detlength/2., zpos + detlength/2.,TrackerThickness);
      geo->AddLayerGeom(layer, mygeom);
      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
      steppingAction_ = new PHG4CylinderSteppingAction(detector_);
      steppingAction_->set_zmin(zpos-detlength/2.);
      steppingAction_->set_zmax(zpos + detlength/2.);
    }
  if (blackhole && !active)
    {
      steppingAction_ = new PHG4CylinderSteppingAction(detector_);
      steppingAction_->set_zmin(zpos-detlength/2.);
      steppingAction_->set_zmax(zpos + detlength/2.);
    }
  if (steppingAction_)
    {
      trackingAction_ = new PHG4FlushStepTrackingAction(steppingAction_);
    }
  return 0;

}

//_______________________________________________________________________
int PHG4CylinderSubsystem::process_event( PHCompositeNode* topNode )
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
PHG4Detector* PHG4CylinderSubsystem::GetDetector( void ) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4CylinderSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}

PHG4TrackingAction*
PHG4CylinderSubsystem::GetTrackingAction( void ) const
{
  return trackingAction_; 
}
