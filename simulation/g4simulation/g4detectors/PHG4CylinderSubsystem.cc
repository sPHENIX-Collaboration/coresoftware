#include "PHG4CylinderSubsystem.h"
#include "PHG4CylinderDetector.h"
#include "PHG4CylinderGeomv1.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4FlushStepTrackingAction.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4CylinderSubsystem::PHG4CylinderSubsystem( const std::string &na, const int lyr):
  PHG4DetectorSubsystem(na,lyr),
  detector_( NULL ),
  steppingAction_( NULL ),
  trackingAction_(NULL),
  eventAction_(NULL)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4CylinderSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  // create hit list only for active layers
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));
  // create detector
  detector_ = new PHG4CylinderDetector(topNode, GetParams(), Name(), GetLayer());
  if (GetParams()->get_int_param("lengthviarapidity"))
    {
      GetParams()->set_double_param("length",PHG4Utils::GetLengthForRapidityCoverage( GetParams()->get_double_param("radius") + GetParams()->get_double_param("thickness"))*2);
    }
  G4double detlength = GetParams()->get_double_param("length");
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
    {
      ostringstream nodename;
      ostringstream geonode;
      if (SuperDetector() != "NONE")
        {
          nodename <<  "G4HIT_" << SuperDetector();
          geonode << "CYLINDERGEOM_" << SuperDetector();
        }
      else
        {
          nodename <<  "G4HIT_" << Name();
          geonode << "CYLINDERGEOM_" << Name();
        }
      PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str() );
      if ( !cylinder_hits )
        {
          dstNode->addNode( new PHIODataNode<PHObject>( cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));
        }
      cylinder_hits->AddLayer(GetLayer());
      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonode.str().c_str());
      if (!geo)
        {
          geo = new PHG4CylinderGeomContainer();
          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
          runNode->addNode(newNode);
        }
      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv1(GetParams()->get_double_param("radius"), GetParams()->get_double_param("place_z")-detlength/2., GetParams()->get_double_param("place_z") + detlength/2.,GetParams()->get_double_param("thickness"));
      geo->AddLayerGeom(GetLayer(), mygeom);
      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
      steppingAction_ = new PHG4CylinderSteppingAction(detector_, GetParams());
    }
  if (GetParams()->get_int_param("blackhole"))
    {
      steppingAction_ = new PHG4CylinderSteppingAction(detector_, GetParams());
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

void
PHG4CylinderSubsystem::SetDefaultParameters()
{
  set_default_double_param("length",100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("radius", 100);
  set_default_double_param("thickness",100);

  set_default_int_param("lengthviarapidity",1);

  set_default_string_param("material", "G4_Air");
}

void
PHG4CylinderSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  GetParams()->Print();
  return;
}
