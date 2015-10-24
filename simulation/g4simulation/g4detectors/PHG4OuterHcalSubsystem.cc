#include "PHG4OuterHcalSubsystem.h"
#include "PHG4OuterHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4OuterHcalSteppingAction.h"
#include "PHG4OuterHcalParameters.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4OuterHcalSubsystem::PHG4OuterHcalSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( NULL ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  detector_type(name),
  superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  params = new PHG4OuterHcalParameters();
}

//_______________________________________________________________________
int PHG4OuterHcalSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4OuterHcalDetector(topNode, params, Name());
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  set<string> nodes;
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
      nodes.insert(nodename.str());
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
          nodes.insert(nodename.str());
	}
      BOOST_FOREACH(string node, nodes)
	{
	  PHG4HitContainer* g4_hits =  findNode::getClass<PHG4HitContainer>( topNode , node.c_str());
	  if ( !g4_hits )
	    {
	      dstNode->addNode( new PHIODataNode<PHObject>( g4_hits = new PHG4HitContainer(), node.c_str(), "PHObject" ));
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
      steppingAction_ = new PHG4OuterHcalSteppingAction(detector_, params);
    }
  if (params->blackhole && !params->active)
    {
      steppingAction_ = new PHG4OuterHcalSteppingAction(detector_, params);
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4OuterHcalSubsystem::process_event( PHCompositeNode * topNode )
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
PHG4OuterHcalSubsystem::Print(const string &what) const
{
  cout << "Outer Hcal Parameters: " << endl;
  params->print();
  detector_->Print(what);
  return;
}


//_______________________________________________________________________
PHG4Detector* PHG4OuterHcalSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4OuterHcalSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

void
PHG4OuterHcalSubsystem::SetTiltAngle(const double tilt)
{
  params->tilt_angle = tilt * deg;
  params->ncross = 0;
}

double
PHG4OuterHcalSubsystem::GetTiltAngle() const
{
  return params->tilt_angle/deg;
}

void
PHG4OuterHcalSubsystem::SetPlaceZ(const G4double dbl)
{
  params->place_in_z = dbl;
}

void
PHG4OuterHcalSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
{
  params->place_in_x = place_x * cm;
  params->place_in_y = place_y * cm;
  params->place_in_z = place_z * cm;
}

void
PHG4OuterHcalSubsystem::SetXRot(const G4double dbl)
{
  params->x_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetYRot(const G4double dbl)
{
  params->y_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetZRot(const G4double dbl)
{
  params->z_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetMaterial(const std::string &mat)
{
  params->material = mat;
}

void
PHG4OuterHcalSubsystem::SetActive(const int i)
{
  params->active = i;
}

void
PHG4OuterHcalSubsystem::SetAbsorberActive(const int i)
{
  params->absorberactive = i;
}

void
PHG4OuterHcalSubsystem::BlackHole(const int i)
{
  params->blackhole = i;
}

void
PHG4OuterHcalSubsystem::SetInnerRadius(const double inner)
{
  params->inner_radius = inner * cm;
}

double
PHG4OuterHcalSubsystem::GetInnerRadius() const
{
  return params->inner_radius/cm;
}

void
PHG4OuterHcalSubsystem::SetOuterRadius(const double outer)
{
  params->outer_radius = outer * cm;
}

double
PHG4OuterHcalSubsystem::GetOuterRadius() const
{
  return params->outer_radius/cm;
}

void
PHG4OuterHcalSubsystem::SetLength(const double len)
{
  params->size_z = len * cm;
}

void
PHG4OuterHcalSubsystem::SetGapWidth(const double gap)
{
  params->scinti_gap = gap *cm;
}

void
PHG4OuterHcalSubsystem::SetNumScintiPlates(const int nplates)
{
  params->n_scinti_plates = nplates;
}

void
PHG4OuterHcalSubsystem::SetNumScintiTiles(const int ntiles)
{
  params->n_scinti_tiles = ntiles;
}

void
PHG4OuterHcalSubsystem::SetScintiThickness(const double thick)
{
  params->scinti_tile_thickness = thick * cm;
}

void
PHG4OuterHcalSubsystem::SetScintiGap(const double scgap)
{
  params->scinti_gap_neighbor = scgap * cm;
}

void
PHG4OuterHcalSubsystem::SetTiltViaNcross(const int ncross)
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
  params->ncross = ncross;
}

void
PHG4OuterHcalSubsystem::SetStepLimits(const double slim) 
{
  params->steplimits = slim*cm;
}

void
PHG4OuterHcalSubsystem::SetLightCorrection(const float inner_radius, const float inner_corr,
			  const float outer_radius, const float outer_corr) 
{
  params->light_balance = true;
  params->light_balance_inner_radius = inner_radius;
  params->light_balance_inner_corr = inner_corr;
  params->light_balance_outer_radius = outer_radius;
  params->light_balance_outer_corr = outer_corr;
}

void
PHG4OuterHcalSubsystem:: SetLightScintModel(const bool b)
{
  params->light_scint_model = b;
}
