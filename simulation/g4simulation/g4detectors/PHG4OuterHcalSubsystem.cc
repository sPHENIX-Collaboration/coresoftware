#include "PHG4OuterHcalSubsystem.h"
#include "PHG4OuterHcalDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4OuterHcalSteppingAction.h"
#include "PHG4OuterHcalParameters.h"

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
  detector_type(name),
  superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  paramsold = new PHG4OuterHcalParameters();
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
  detector_ = new PHG4OuterHcalDetector(topNode, paramsold, params, Name());
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(overlapcheck);
  set<string> nodes;
  if (paramsold->active)
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
      if (paramsold->absorberactive)
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
      steppingAction_->EnableFieldChecker(enable_field_checker);
    }
  if (paramsold->blackhole && !paramsold->active)
    {
      steppingAction_ = new PHG4OuterHcalSteppingAction(detector_, params);
      steppingAction_->EnableFieldChecker(enable_field_checker);
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
PHG4OuterHcalSubsystem::SetTiltAngle(const double tilt)
{
  paramsold->tilt_angle = tilt * deg;
  paramsold->ncross = 0;
}

double
PHG4OuterHcalSubsystem::GetTiltAngle() const
{
  return paramsold->tilt_angle/deg;
}

void
PHG4OuterHcalSubsystem::SetPlaceZ(const G4double dbl)
{
  paramsold->place_in_z = dbl;
}

void
PHG4OuterHcalSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
{
  paramsold->place_in_x = place_x * cm;
  paramsold->place_in_y = place_y * cm;
  paramsold->place_in_z = place_z * cm;
}

void
PHG4OuterHcalSubsystem::SetXRot(const G4double dbl)
{
  paramsold->x_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetYRot(const G4double dbl)
{
  paramsold->y_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetZRot(const G4double dbl)
{
  paramsold->z_rot = dbl * deg;
}

void
PHG4OuterHcalSubsystem::SetMaterial(const std::string &mat)
{
  paramsold->material = mat;
}

void
PHG4OuterHcalSubsystem::SetActive(const int i)
{
  paramsold->active = i;
}

void
PHG4OuterHcalSubsystem::SetAbsorberActive(const int i)
{
  paramsold->absorberactive = i;
}

void
PHG4OuterHcalSubsystem::BlackHole(const int i)
{
  paramsold->blackhole = i;
}

void
PHG4OuterHcalSubsystem::SetInnerRadius(const double inner)
{
  paramsold->inner_radius = inner * cm;
}

double
PHG4OuterHcalSubsystem::GetInnerRadius() const
{
  return paramsold->inner_radius/cm;
}

void
PHG4OuterHcalSubsystem::SetOuterRadius(const double outer)
{
  paramsold->outer_radius = outer * cm;
}

double
PHG4OuterHcalSubsystem::GetOuterRadius() const
{
  return paramsold->outer_radius/cm;
}

void
PHG4OuterHcalSubsystem::SetLength(const double len)
{
  paramsold->size_z = len * cm;
}

void
PHG4OuterHcalSubsystem::SetGapWidth(const double gap)
{
  paramsold->scinti_gap = gap *cm;
}

void
PHG4OuterHcalSubsystem::SetNumScintiPlates(const int nplates)
{
  paramsold->n_scinti_plates = nplates;
}

void
PHG4OuterHcalSubsystem::SetNumScintiTiles(const int ntiles)
{
  paramsold->n_scinti_tiles = ntiles;
}

void
PHG4OuterHcalSubsystem::SetScintiThickness(const double thick)
{
  paramsold->scinti_tile_thickness = thick * cm;
}

void
PHG4OuterHcalSubsystem::SetScintiGap(const double scgap)
{
  paramsold->scinti_gap_neighbor = scgap * cm;
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
  paramsold->ncross = ncross;
}

void
PHG4OuterHcalSubsystem::SetStepLimits(const double slim) 
{
  paramsold->steplimits = slim*cm;
}

void
PHG4OuterHcalSubsystem::SetLightCorrection(const float inner_radius, const float inner_corr,
			  const float outer_radius, const float outer_corr) 
{
  paramsold->light_balance = true;
  paramsold->light_balance_inner_radius = inner_radius;
  paramsold->light_balance_inner_corr = inner_corr;
  paramsold->light_balance_outer_radius = outer_radius;
  paramsold->light_balance_outer_corr = outer_corr;
}

void
PHG4OuterHcalSubsystem:: SetLightScintModel(const bool b)
{
  paramsold->light_scint_model = b;
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

