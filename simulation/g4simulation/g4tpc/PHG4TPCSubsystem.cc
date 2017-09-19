#include "PHG4TPCSubsystem.h"
#include "PHG4TPCDetector.h"
#include "PHG4TPCSteppingAction.h"

#include <g4detectors/PHG4Parameters.h>

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <set>
#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4TPCSubsystem::PHG4TPCSubsystem( const std::string &name, const int lyr ):
  PHG4DetectorSubsystem( name, lyr ),
  detector_(nullptr),
  steppingAction_( nullptr )
{
  InitializeParameters();
}

//_______________________________________________________________________
int 
PHG4TPCSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4TPCDetector(topNode, GetParams(), Name());
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());
  set<string> nodes;
  if (GetParams()->get_int_param("active"))
    {
      PHNodeIterator dstIter( dstNode );
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode",SuperDetector()));
      if (! DetNode)
	{
          DetNode = new PHCompositeNode(SuperDetector());
          dstNode->addNode(DetNode);
        }

      ostringstream nodename;
      if (SuperDetector() != "NONE")
	{
	  nodename <<  "G4HIT_" << SuperDetector();
	}
      else
	{
	  nodename <<  "G4HIT_" << Name();
	}
      nodes.insert(nodename.str());
      if (GetParams()->get_int_param("absorberactive"))
	{
	  nodename.str("");
	  if (SuperDetector() != "NONE")
	    {
	      nodename <<  "G4HIT_ABSORBER_" << SuperDetector();
	    }
	  else
	    {
	      nodename <<  "G4HIT_ABSORBER_" << Name();
	    }
          nodes.insert(nodename.str());
	}
      BOOST_FOREACH(string node, nodes)
	{
	  PHG4HitContainer* g4_hits =  findNode::getClass<PHG4HitContainer>( topNode , node.c_str());
	  if ( !g4_hits )
	    {
	      g4_hits = new PHG4HitContainer(node);
              DetNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));
	    }
	}

      // create stepping action
      steppingAction_ = new PHG4TPCSteppingAction(detector_, GetParams());
    }
  else
    {
      // if this is a black hole it does not have to be active
      if (GetParams()->get_int_param("blackhole"))
	{
	  steppingAction_ = new PHG4TPCSteppingAction(detector_, GetParams());
	}
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4TPCSubsystem::process_event( PHCompositeNode * topNode )
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
PHG4TPCSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  GetParams()->Print();
  if (detector_)
    {
      detector_->Print(what);
    }
  if (steppingAction_)
    {
      steppingAction_->Print(what);
    }

  return;
}

//_______________________________________________________________________
PHG4Detector* PHG4TPCSubsystem::GetDetector( void ) const
{
  return detector_;
}

void
PHG4TPCSubsystem::SetDefaultParameters()
{
  set_default_double_param("gas_inner_radius",21.);
  set_default_double_param("gas_outer_radius",77.);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("thickness_cu",0.005);
  set_default_double_param("thickness_fr4",0.005);
  set_default_double_param("thickness_honeycomb",1.);
  set_default_double_param("thickness_kapton",0.15);
  set_default_double_param("thickness_pcb",0.005);
  set_default_double_param("tpc_length",211.);

  set_default_double_param("steplimits", NAN);

  set_default_string_param("tpc_gas", "sPHENIX_TPC_Gas");
}


