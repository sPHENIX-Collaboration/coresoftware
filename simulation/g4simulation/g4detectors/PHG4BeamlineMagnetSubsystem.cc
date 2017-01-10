#include "PHG4BeamlineMagnetSubsystem.h"
#include "PHG4BeamlineMagnetDetector.h"
#include "PHG4CylinderGeomv1.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4BeamlineMagnetSubsystem::PHG4BeamlineMagnetSubsystem( const std::string &na, const int lyr):
  PHG4DetectorSubsystem(na,lyr),
  detector_( NULL ),
  steppingAction_( NULL ),
  eventAction_(NULL)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4BeamlineMagnetSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  // create hit list only for active layers
  PHNodeIterator iter( topNode );
  //  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));
  if (GetParams()->get_int_param("lengthviarapidity"))
    {
      GetParams()->set_double_param("length",PHG4Utils::GetLengthForRapidityCoverage( GetParams()->get_double_param("radius") + GetParams()->get_double_param("thickness"))*2);
    }
  // create detector
  detector_ = new PHG4BeamlineMagnetDetector(topNode, GetParams(), Name(), GetLayer());
//  G4double detlength = GetParams()->get_double_param("length");
//  detector_->SuperDetector(SuperDetector());
//  detector_->OverlapCheck(CheckOverlap());
//  if (GetParams()->get_int_param("active"))
//    {
//      ostringstream nodename;
//      ostringstream geonode;
//      if (SuperDetector() != "NONE")
//        {
//          nodename <<  "G4HIT_" << SuperDetector();
//          geonode << "CYLINDERGEOM_" << SuperDetector();
//        }
//      else
//        {
//          nodename <<  "G4HIT_" << Name();
//          geonode << "CYLINDERGEOM_" << Name();
//        }
//      PHG4HitContainer* cylinder_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str() );
//      if ( !cylinder_hits )
//        {
//          dstNode->addNode( new PHIODataNode<PHObject>( cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));
//        }
//      cylinder_hits->AddLayer(GetLayer());
//      PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonode.str().c_str());
//      if (!geo)
//        {
//          geo = new PHG4CylinderGeomContainer();
//          PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
//          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo, geonode.str().c_str(), "PHObject");
//          runNode->addNode(newNode);
//        }
//      PHG4CylinderGeom *mygeom = new PHG4CylinderGeomv1(GetParams()->get_double_param("radius"), GetParams()->get_double_param("place_z")-detlength/2., GetParams()->get_double_param("place_z") + detlength/2.,GetParams()->get_double_param("thickness"));
//      geo->AddLayerGeom(GetLayer(), mygeom);
//      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
//      steppingAction_ = new PHG4CylinderSteppingAction(detector_, GetParams());
//    }
//  if (GetParams()->get_int_param("blackhole"))
//    {
//      steppingAction_ = new PHG4CylinderSteppingAction(detector_, GetParams());
//    }
  return 0;

}

//_______________________________________________________________________
int PHG4BeamlineMagnetSubsystem::process_event( PHCompositeNode* topNode )
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
PHG4Detector* PHG4BeamlineMagnetSubsystem::GetDetector( void ) const
{
  return detector_;
}

void
PHG4BeamlineMagnetSubsystem::SetDefaultParameters()
{
  set_default_string_param("magtype", "");

  set_default_double_param("field_y", NAN);
  set_default_double_param("fieldgradient", NAN);

  set_default_double_param("length",100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("radius", 100);
  set_default_double_param("thickness",100);
  set_default_double_param("tmin",NAN);
  set_default_double_param("tmax",NAN);

  set_default_int_param("lengthviarapidity",1);

  set_default_string_param("material", "G4_Galactic");
}

void
PHG4BeamlineMagnetSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  if (! BeginRunExecuted())
    {
      cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
      cout << "To do so either run one or more events or on the command line execute: " << endl;
      cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
      cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
      cout << "g4->InitRun(se->topNode());" << endl;
      cout << "PHG4BeamlineMagnetSubsystem *cyl = (PHG4BeamlineMagnetSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
      cout << "cyl->Print()" << endl;
      return;
    }
  GetParams()->Print();
  return;
}
