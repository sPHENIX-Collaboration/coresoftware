/*===============================================================*
 *                        March 2nd 2017                         *
 *         mRICH Subsystem created by Cheuk-Ping Wong @GSU       *
 *===============================================================*/
#include "PHG4mRICHSubsystem.h"
#include "PHG4mRICHDetector.h"
#include "PHG4Parameters.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4mRICHSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>
#include <sstream>
#include <boost/foreach.hpp>
//#include <string>

using namespace std;

//_______________________________________________________________________
PHG4mRICHSubsystem::PHG4mRICHSubsystem( const std::string &name, const int lyr):
  PHG4DetectorSubsystem( name, lyr ),
  _detector( NULL ),
  _detectorName(name),
  _steppingAction(NULL),
  _eventAction(NULL)
{
  InitializeParameters();
}

//_______________________________________________________________________
int
PHG4mRICHSubsystem::InitSubsystem( PHCompositeNode* topNode )
{
  // kludge until the phg4parameters are sorted out (adding layers)
  GetParams()->set_name(Name());
  GetParams()->set_int_param("active",1);
  GetParams()->set_int_param("blackhole",1);
  GetParams()->set_int_param("absorberactive",1);
  GetParams()->set_string_param("superdetector","NONE");
  GetParams()->set_string_param("detectorname","mRICH");
  return 0;
}

//_______________________________________________________________________
int PHG4mRICHSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  //---------------------------------
  // create detector
  //---------------------------------
  _detector = new PHG4mRICHDetector(topNode, GetParams(), Name(), GetLayer());
  _detector->SuperDetector(SuperDetector());
  _detector->OverlapCheck(CheckOverlap());
  
  //---------------------------------
  // create hit node and stepping action
  //---------------------------------
  if (GetParams()->get_int_param("active")) {
    set<string> nodes;
    
    // create hit output node
    ostringstream nodename;
    nodename <<  "G4HIT_" << GetParams()->get_string_param("detectorname");
    
    PHG4HitContainer* mRICH_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
    if (!mRICH_hits) {
      mRICH_hits = new PHG4HitContainer(nodename.str());
      PHIODataNode<PHObject> *hitNode = new PHIODataNode<PHObject>(mRICH_hits, nodename.str().c_str(), "PHObject");
      dstNode->addNode(hitNode);
      nodes.insert(nodename.str());
    }
    
    ostringstream absnodename;
    absnodename << "G4HIT_ABSORBER_" << GetParams()->get_string_param("detectorname");
    
    PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
    if (!absorber_hits) {
      absorber_hits = new PHG4HitContainer(absnodename.str());
      PHIODataNode<PHObject> *abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
      dstNode->addNode(abshitNode);
      nodes.insert(nodename.str());
    }
    
    // create stepping action
    _steppingAction = new PHG4mRICHSteppingAction(_detector,GetParams());
    
    // event actions
    BOOST_FOREACH(string node, nodes) {
      if (!_eventAction) _eventAction = new PHG4EventActionClearZeroEdep(topNode, node);
      else {
	PHG4EventActionClearZeroEdep *evtact =
	  dynamic_cast<PHG4EventActionClearZeroEdep *>(_eventAction);
	evtact->AddNode(node);
      }
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4mRICHSubsystem::process_event( PHCompositeNode* topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if(_steppingAction)  _steppingAction->SetInterfacePointers( topNode );
  return 0;
}


//_______________________________________________________________________
PHG4Detector*
PHG4mRICHSubsystem::GetDetector( void ) const
{
  return _detector;
}
//_______________________________________________________________________
void PHG4mRICHSubsystem::SetDefaultParameters()
{
  set_default_int_param("single_mRICH", 1);    //1 for single mRICH
                                               //0 for mRICH wall

  set_default_int_param("use_g4steps",0);        //for stepping function
}
