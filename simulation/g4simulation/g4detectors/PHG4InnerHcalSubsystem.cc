#include "PHG4InnerHcalSubsystem.h"
#include "PHG4InnerHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4InnerHcalSteppingAction.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <set>
#include <sstream>

using namespace std;

map<const string,double> default_double;
map<const string, int> default_int;
map<const string, string> default_string;


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
}

int 
PHG4InnerHcalSubsystem::Init(PHCompositeNode* topNode)
{
  params = new PHG4Parameters(superdetector);
  SetDefaultParameters();
  return 0;
}

//_______________________________________________________________________
int 
PHG4InnerHcalSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  string g4geonodename = "G4GEO_" + superdetector;
  string paramnodename = "G4GEOPARAM_" + superdetector;

  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode,paramnodename);
  if (nodeparams)
    {
      params->FillFrom(nodeparams);
    }
  else
    {
      params->SaveToNodeTree(parNode,paramnodename);
    }
  UpdateParametersWithMacro();
  parNode->addNode(new PHDataNode<PHG4Parameters>(params,g4geonodename));
  // create detector
  detector_ = new PHG4InnerHcalDetector(topNode, params, Name());
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  set<string> nodes;
  if (params->get_int_param("active"))
    {
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",superdetector));
      if (! DetNode)
	{
          DetNode = new PHCompositeNode(superdetector);
          dstNode->addNode(DetNode);
        }

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
      if (params->get_int_param("absorberactive"))
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
	      g4_hits = new PHG4HitContainer();
	      //	      DetNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));
	      dstNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));

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
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, params);

    }
  else
    {
      // if this is a black hole it does not have to be active
      if (params->get_int_param("blackhole"))
	{
	  steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, params);
	}
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


void
PHG4InnerHcalSubsystem::Print(const string &what) const
{
  cout << "Inner Hcal Parameters: " << endl;
  params->print();
  if (detector_)
    {
      detector_->Print(what);
    }
  return;
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
PHG4InnerHcalSubsystem::SetActive(const int i)
{
  iparams["active"] = i;
}

void
PHG4InnerHcalSubsystem::SetAbsorberActive(const int i)
{
  iparams["absorberactive"] = i;
}

void
PHG4InnerHcalSubsystem::BlackHole(const int i)
{
  iparams["blackhole"] = i;
}

void
PHG4InnerHcalSubsystem::set_double_param(const std::string &name, const double dval)
{
  if (default_double.find(name) == default_double.end())
    {
      cout << "double parameter " << name << " not implemented" << endl;
      cout << "implemented double parameters are:" << endl;
      for (map<const string, double>::const_iterator iter = default_double.begin(); iter != default_double.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  dparams[name] = dval;
}

void
PHG4InnerHcalSubsystem::set_int_param(const std::string &name, const int ival)
{
  if (default_int.find(name) == default_int.end())
    {
      cout << "integer parameter " << name << " not implemented" << endl;
      cout << "implemented integer parameters are:" << endl;
      for (map<const string, int>::const_iterator iter = default_int.begin(); iter != default_int.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  iparams[name] = ival;
}

void
PHG4InnerHcalSubsystem::set_string_param(const std::string &name, const string &sval)
{
  if (default_string.find(name) == default_string.end())
    {
      cout << "string parameter " << name << " not implemented" << endl;
      cout << "implemented string parameters are:" << endl;
      for (map<const string, string>::const_iterator iter = default_string.begin(); iter != default_string.end(); ++iter)
	{
	  cout << iter->first << endl;
	}
      return;
    }
  cparams[name] = sval;
}

void
PHG4InnerHcalSubsystem::SetDefaultParameters()
{

  default_double["inner_radius"] = 116.;
  default_double["light_balance_inner_corr"] = NAN;
  default_double["light_balance_inner_radius"] = NAN;
  default_double["light_balance_outer_corr"] = NAN;
  default_double["light_balance_outer_radius"] = NAN;
  default_double["outer_radius"] = 136.;
  default_double["place_x"] = 0.;
  default_double["place_y"] = 0.;
  default_double["place_z"] = 0.;
  default_double["rot_x"] = 0.;
  default_double["rot_y"] = 0.;
  default_double["rot_z"] = 0.;
  default_double["scinti_eta_coverage"] = 1.1;
  default_double["scinti_gap"] = 0.85;
  default_double["scinti_gap_neighbor"] = 0.1;
  default_double["scinti_tile_thickness"] = 0.7;
  default_double["size_z"] = 175.94 * 2;
  default_double["steplimits"] = NAN;
  default_double["tilt_angle"] = NAN; // default is 4 crossinge
  default_double["timing_cut"] = 100.;

  default_int["absorberactive"] = 0;
  default_int["absorbertruth"] = 0;
  default_int["active"] = 0;
  default_int["blackhole"] = 0;
  default_int["light_scint_model"] = 1;
  default_int["ncross"] = 4;
  default_int["n_scinti_plates"] = 5*64;
  default_int["n_scinti_tiles"] = 12;

  default_string["material"] = "SS310";
  for (map<const string,double>::const_iterator iter = default_double.begin(); iter != default_double.end(); ++iter)
    {
      params->set_double_param(iter->first,iter->second);
    }
  for (map<const string,int>::const_iterator iter = default_int.begin(); iter != default_int.end(); ++iter)
    {
      params->set_int_param(iter->first,iter->second);
    }
  for (map<const string,string>::const_iterator iter = default_string.begin(); iter != default_string.end(); ++iter)
    {
      params->set_string_param(iter->first,iter->second);
    }

}

void
PHG4InnerHcalSubsystem::UpdateParametersWithMacro()
{
  for (map<const string,double>::const_iterator iter = dparams.begin(); iter != dparams.end(); ++iter)
    {
      params->set_double_param(iter->first,iter->second);
    }
  for (map<const string,int>::const_iterator iter = iparams.begin(); iter != iparams.end(); ++iter)
    {
      params->set_int_param(iter->first,iter->second);
    }
  for (map<const string,string>::const_iterator iter = cparams.begin(); iter != cparams.end(); ++iter)
    {
      params->set_string_param(iter->first,iter->second);
    }
  return;
}

