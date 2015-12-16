#include "PHG4OuterHcalSubsystem.h"
#include "PHG4OuterHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4OuterHcalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

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
  enable_field_checker(0),
  layer(lyr),
  usedb(0),
  filetype(PHG4OuterHcalSubsystem::none),
  detector_type(name),
  superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  params = new PHG4Parameters(Name()); // temporary name till the init is called
  SetDefaultParameters();  
}

void
PHG4OuterHcalSubsystem::SuperDetector(const std::string &name)
{
  superdetector = name;
  Name(name);
  return;
}

int 
PHG4OuterHcalSubsystem::Init(PHCompositeNode* topNode)
{
  params->set_name(superdetector);
  return 0;
}

//_______________________________________________________________________
int PHG4OuterHcalSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  string g4geonodename = "G4GEO_" + superdetector;
  parNode->addNode(new PHDataNode<PHG4Parameters>(params,g4geonodename));

  string paramnodename = "G4GEOPARAM_" + superdetector;
  // ASSUMPTION: if we read from DB and/or file we don't want the stuff from
  // the node tree
  // We leave the defaults intact in case there is no entry for
  // those in the object read from the DB or file
  // Order: read first DB, then calib file if both are enabled
  if (usedb || filetype != PHG4OuterHcalSubsystem::none)
    {
      if (usedb)
	{
          ReadParamsFromDB();
	}
      if (filetype != PHG4OuterHcalSubsystem::none)
	{
	  ReadParamsFromFile(filetype);
	}
    }
  else
    {
      PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode,paramnodename);
      if (nodeparams)
	{
	  params->FillFrom(nodeparams);
	}
    }
  // parameters set in the macro always override whatever is read from
  // the node tree, DB or file
  UpdateParametersWithMacro();
  // save updated persistant copy on node tree
  params->SaveToNodeTree(parNode,paramnodename);
  // create detector
  detector_ = new PHG4OuterHcalDetector(topNode, params, Name());
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
	      DetNode->addNode( new PHIODataNode<PHObject>( g4_hits, node.c_str(), "PHObject" ));
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
      steppingAction_->EnableFieldChecker(enable_field_checker);
    }
  else
    {
      if (params->get_int_param("blackhole"))
	{
	  steppingAction_ = new PHG4OuterHcalSteppingAction(detector_, params);
	  steppingAction_->EnableFieldChecker(enable_field_checker);
	}
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
  if (detector_)
    {
      detector_->Print(what);
    }
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
PHG4OuterHcalSubsystem::SetActive(const int i)
{
  iparams["active"] = i;
}

void
PHG4OuterHcalSubsystem::SetAbsorberActive(const int i)
{
  iparams["absorberactive"] = i;
}

void
PHG4OuterHcalSubsystem::BlackHole(const int i)
{
  iparams["blackhole"] = i;
}

void
PHG4OuterHcalSubsystem::SetLightCorrection(const double inner_radius, const double inner_corr,const double outer_radius, const double outer_corr)
{
  dparams["light_balance_inner_corr"] = inner_corr;
  dparams["light_balance_inner_radius"] = inner_radius;
  dparams["light_balance_outer_corr"] = outer_corr;
  dparams["light_balance_outer_radius"] = outer_radius;
  return;
}

void
PHG4OuterHcalSubsystem::set_double_param(const std::string &name, const double dval)
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
PHG4OuterHcalSubsystem::set_int_param(const std::string &name, const int ival)
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
PHG4OuterHcalSubsystem::set_string_param(const std::string &name, const string &sval)
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
PHG4OuterHcalSubsystem::SetDefaultParameters()
{

  default_double["inner_radius"] = 178.;
  default_double["light_balance_inner_corr"] = NAN;
  default_double["light_balance_inner_radius"] = NAN;
  default_double["light_balance_outer_corr"] = NAN;
  default_double["light_balance_outer_radius"] = NAN;
  default_double["magnet_cutout"] = 12.;
  default_double["outer_radius"] = 260.;
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
  default_double["size_z"] = 304.91 * 2;
  default_double["steplimits"] = NAN;
  default_double["tilt_angle"] = NAN; // default is 4 crossinge

  default_int["absorberactive"] = 0;
  default_int["absorbertruth"] = 0;
  default_int["active"] = 0;
  default_int["blackhole"] = 0;
  default_int["light_scint_model"] = 1;
  default_int["magnet_cutout_first_scinti"] = 8; // tile start at 0, drawing tile starts at 1
  default_int["ncross"] = -4;
  default_int["n_scinti_plates"] = 5*64;
  default_int["n_scinti_tiles"] = 12;

  default_string["material"] = "G4_Fe";

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
PHG4OuterHcalSubsystem::UpdateParametersWithMacro()
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

int
PHG4OuterHcalSubsystem::SaveParamsToDB()
{
  int iret = params->WriteToDB();
  if (iret)
    {
      cout << "problem committing to DB" << endl;
    }
  return iret;
}

int
PHG4OuterHcalSubsystem::ReadParamsFromDB()
{
  int iret = params->ReadFromDB();
  if (iret)
    {
      cout << "problem reading from DB" << endl;
    }
  params->print();
  return iret;
}

int
PHG4OuterHcalSubsystem::SaveParamsToFile(const PHG4OuterHcalSubsystem::FILE_TYPE ftyp)
{
  string extension;
  switch(ftyp)
    {
    case xml:
      extension = "xml";
      break;
    case root:
      extension = "root";
      break;
    default:
      cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
      exit(1);
    }

  int iret = params->WriteToFile(extension,calibfiledir);
  if (iret)
    {
      cout << "problem saving to " << extension << " file " << endl;
    }
  return iret;
}

int
PHG4OuterHcalSubsystem::ReadParamsFromFile(const PHG4OuterHcalSubsystem::FILE_TYPE ftyp)
{
  string extension;
  switch(ftyp)
    {
    case xml:
      extension = "xml";
      break;
    case root:
      extension = "root";
      break;
    default:
      cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
      exit(1);
    }
  int iret = params->ReadFromFile(extension,calibfiledir);
  if (iret)
    {
      cout << "problem saving to " << extension << " file " << endl;
    }
  return iret;
}


double
PHG4OuterHcalSubsystem::get_double_param(const std::string &name) const
{
  return params->get_double_param(name);
}

int
PHG4OuterHcalSubsystem::get_int_param(const std::string &name) const
{
  return params->get_int_param(name);
}

string
PHG4OuterHcalSubsystem::get_string_param(const std::string &name) const
{
  return params->get_string_param(name);
}

