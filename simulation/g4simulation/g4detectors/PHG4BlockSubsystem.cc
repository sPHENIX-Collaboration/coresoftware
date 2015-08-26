#include "PHG4BlockSubsystem.h"
#include "PHG4BlockDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4BlockRegionSteppingAction.h"
#include "PHG4BlockSteppingAction.h"
#include "PHG4BlockGeomv1.h"
#include "PHG4BlockGeomContainer.h"
#include <g4main/PHG4Utils.h>

#include <g4main/PHG4HitContainer.h>
#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4BlockSubsystem::PHG4BlockSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  _detector( 0 ),
  _steppingAction(NULL),
  _eventAction(NULL),
  _center_in_x(0),
  _center_in_y(0),
  _center_in_z(0),
  _rot_in_z(0),
  _material("G4_Galactic"),  // default - almost nothing
  _active(0),
  _layer(lyr),
  _blackhole(0),
  _use_g4_steps(0),
  _use_ionisation_energy(0),
  _detector_type(name),
  _superdetector("NONE")
{

  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  for (int i = 0; i < 3; i++) {
    _dimension[i] = 100.0 * cm;
  }
}

//_______________________________________________________________________
int PHG4BlockSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  _detector = new PHG4BlockDetector(topNode, Name(), _layer);
  _detector->SetSize(_dimension[0], _dimension[1], _dimension[2]);
  _detector->SetCenter(_center_in_x, _center_in_y, _center_in_z);
  _detector->SetZRot(_rot_in_z);
  _detector->SetMaterial(_material);
  _detector->SetActive(_active);
  _detector->BlackHole(_blackhole);
  _detector->SuperDetector(_superdetector);
  _detector->OverlapCheck(overlapcheck);
  if(_active) {

    ostringstream nodename;
    ostringstream geonode;
    if(_superdetector != "NONE") {
      nodename <<  "G4HIT_" << _superdetector;
      geonode << "BLOCKGEOM_" << _superdetector;
    } else {
      nodename <<  "G4HIT_" << _detector_type << "_" << _layer;
      geonode << "BLOCKGEOM_" << _detector_type << "_" << _layer;
    }

    // create hit list
    PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
    if( !block_hits ){
      dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));
    }

    block_hits->AddLayer(_layer);
    PHG4BlockGeomContainer *geocont = findNode::getClass<PHG4BlockGeomContainer>(topNode,
                                                                                 geonode.str().c_str());
    if(!geocont) {
      geocont = new PHG4BlockGeomContainer();
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geocont, geonode.str().c_str(), "PHObject");
      runNode->addNode(newNode);
    }

    PHG4BlockGeom *geom = new PHG4BlockGeomv1(_layer,
                                              _dimension[0], _dimension[1], _dimension[2],
                                              _center_in_x, _center_in_y, _center_in_z, _rot_in_z);
    geocont->AddLayerGeom(_layer, geom);

    _steppingAction = new PHG4BlockSteppingAction(_detector);
    _steppingAction->UseG4Steps(_use_g4_steps);
    _steppingAction->UseIonizationEnergy(_use_ionisation_energy);
    _eventAction = new PHG4EventActionClearZeroEdep(topNode, nodename.str());

  } else if(_blackhole) {
    _steppingAction = new PHG4BlockSteppingAction(_detector);
  }

  return 0;
}

//_______________________________________________________________________
int PHG4BlockSubsystem::process_event( PHCompositeNode* topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if(_steppingAction) {
    _steppingAction->SetInterfacePointers( topNode );
  }
  return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4BlockSubsystem::GetDetector( void ) const
{
  return _detector;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4BlockSubsystem::GetSteppingAction( void ) const
{
  return _steppingAction;
}

//_______________________________________________________________________
void PHG4BlockSubsystem::UseG4Steps(const int i)
{
  _use_g4_steps = i;
  if(_steppingAction)
  {
    _steppingAction->UseG4Steps(i);
  }
  return;
}

//_______________________________________________________________________
void PHG4BlockSubsystem::UseIonizationEnergy(const int i)
{
  _use_ionisation_energy = i;
  if(_steppingAction)
  {
    _steppingAction->UseIonizationEnergy(i);
  }
  return;
}
