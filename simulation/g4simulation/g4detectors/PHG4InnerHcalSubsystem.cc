#include "PHG4InnerHcalSubsystem.h"
#include "PHG4InnerHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4InnerHcalSteppingAction.h"
#include "PHG4InnerHcalParameters.h"

#include <g4main/PHG4HitContainer.h>
#include <fun4all/getClass.h>

#include <Geant4/globals.hh>

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
int PHG4InnerHcalSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4InnerHcalDetector(topNode, params, Name());
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  if (params->active)
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
      // create hit list
      PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));

	}
      if (params->absorberactive)
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
	}
      block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( new PHG4HitContainer(), nodename.str().c_str(), "PHObject" ));

	}
      // create stepping action
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_);

      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
    }
  if (params->blackhole && !params->active)
    {
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_);
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

void
PHG4InnerHcalSubsystem::SetTiltAngle(const double tilt)
{
  params->tilt_angle = tilt * deg;
}

void
PHG4InnerHcalSubsystem::SetPlaceZ(const G4double dbl)
{
  params->place_in_z = dbl;
}

void
PHG4InnerHcalSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
{
  params->place_in_x = place_x * cm;
  params->place_in_y = place_y * cm;
  params->place_in_z = place_z * cm;
}

void
PHG4InnerHcalSubsystem::SetXRot(const G4double dbl)
{
  params->x_rot = dbl * deg;
}

void
PHG4InnerHcalSubsystem::SetYRot(const G4double dbl)
{
  params->y_rot = dbl * deg;
}

void
PHG4InnerHcalSubsystem::SetZRot(const G4double dbl)
{
  params->z_rot = dbl * deg;
}

void
PHG4InnerHcalSubsystem::SetMaterial(const std::string &mat)
{
  params->material = mat;
}

void
PHG4InnerHcalSubsystem::SetActive(const int i)
{
  params->active = i;
}

void
PHG4InnerHcalSubsystem::SetAbsorberActive(const int i)
{
  params->absorberactive = i;
}

void
PHG4InnerHcalSubsystem::BlackHole(const int i)
{
  params->blackhole = i;
}

void
PHG4InnerHcalSubsystem::SetInnerRadius(const double inner)
{
  params->inner_radius = inner *cm;
}

double
PHG4InnerHcalSubsystem::GetInnerRadius() const
{
  return params->inner_radius/cm;
}

void
PHG4InnerHcalSubsystem::SetOuterRadius(const double outer)
{
  params->outer_radius = outer * cm;
}

double
PHG4InnerHcalSubsystem::GetOuterRadius() const
{
  return params->outer_radius/cm;
}

void
PHG4InnerHcalSubsystem::SetLength(const double len)
{
  params->size_z = len * cm;
}
void
PHG4InnerHcalSubsystem::SetGapWidth(const double gap)
{
  params->scinti_gap = gap *cm;
}

void
PHG4InnerHcalSubsystem::SetNumScintiPlates(const int nplates)
{
  params->n_scinti_plates = nplates;
}

void
PHG4InnerHcalSubsystem::SetNumScintiTiles(const int ntiles)
{
  params->n_scinti_tiles = ntiles;
}

void
PHG4InnerHcalSubsystem::SetScintiThickness(const double thick)
{
  params->scinti_tile_thickness = thick * cm;
}

void
PHG4InnerHcalSubsystem::SetScintiGap(const double scgap)
{
  params->scinti_gap_neighbor = scgap * cm;
}
