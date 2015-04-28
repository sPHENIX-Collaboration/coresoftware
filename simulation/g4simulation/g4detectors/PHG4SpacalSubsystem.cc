// $$Id: PHG4SpacalSubsystem.cc,v 1.2 2014/08/12 03:49:12 jinhuang Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#include "PHG4SpacalSubsystem.h"
#include "PHG4SpacalDetector.h"
#include "PHG4ProjSpacalDetector.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4SpacalSteppingAction.h"
#include "PHG4CylinderEventAction.h"
#include <g4main/PHG4Utils.h>

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4SpacalSubsystem::PHG4SpacalSubsystem( const std::string &na, const int lyr):
  detector_( NULL ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  active(0),
  absorberactive(0),
  layer(lyr),
  lengthViaRapidityCoverage(true),
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
int PHG4SpacalSubsystem::InitRun( PHCompositeNode* topNode )
{
  // create hit list only for active layers
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));
  // create detector
  _geom.set_layer( layer );
  if (lengthViaRapidityCoverage)
    {
       const double half_length =  PHG4Utils::GetLengthForRapidityCoverage(_geom.get_max_radius());
       _geom.set_zmin(-half_length);
       _geom.set_zmax(+half_length);
    }

  switch (_geom.get_config())
    {
  case PHG4CylinderGeom_Spacalv1::kNonProjective:
    cout << "PHG4SpacalSubsystem::InitRun - use PHG4SpacalDetector" << endl;
    detector_ = new PHG4SpacalDetector(topNode, Name(),
        dynamic_cast<PHG4SpacalDetector::SpacalGeom_t *>(&_geom), layer);
    break;

  case PHG4CylinderGeom_Spacalv1::kProjective_PolarTaper:
    cout << "PHG4SpacalSubsystem::InitRun - use PHG4ProjSpacalDetector" << endl;
    detector_ = new PHG4ProjSpacalDetector(topNode, Name(),
        dynamic_cast<PHG4ProjSpacalDetector::SpacalGeom_t *>(&_geom), layer);
    break;

  default:
    cout << "PHG4SpacalSubsystem::InitRun - use PHG4SpacalDetector" << endl;
    exit(1);
    break;

    }

  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);

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
      steppingAction_ = new PHG4SpacalSteppingAction(detector_);
    }


  return 0;

}

//_______________________________________________________________________
int PHG4SpacalSubsystem::process_event( PHCompositeNode* topNode )
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
PHG4Detector* PHG4SpacalSubsystem::GetDetector( void ) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4SpacalSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}


void
PHG4SpacalSubsystem::Print(const std::string &what) const
{
  detector_->Print(what);
  return;
}

