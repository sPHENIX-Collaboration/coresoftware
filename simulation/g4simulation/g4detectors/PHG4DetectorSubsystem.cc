#include "PHG4DetectorSubsystem.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <iostream>
#include <sstream>

using namespace std;

PHG4DetectorSubsystem::PHG4DetectorSubsystem(const std::string &name, const int lyr): 
  PHG4Subsystem(name),
  params(new PHParameters(Name())),
  paramscontainer(NULL),
  savetopNode(NULL),
  overlapcheck(false),
  layer(lyr),
  usedb(0),
  beginrunexecuted(0),
  filetype(PHG4DetectorSubsystem::none),
  superdetector("NONE"),
  calibfiledir("./")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str());
}

int 
PHG4DetectorSubsystem::Init(PHCompositeNode* topNode)
{
  savetopNode = topNode;
  params->set_name(Name());
  int iret = InitSubsystem(topNode);
  return iret;
}

int 
PHG4DetectorSubsystem::InitRun( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));

  string g4geonodename = "G4GEO_";
  string paramnodename = "G4GEOPARAM_";
  string calibdetname;
  int isSuperDetector = 0;
  if (superdetector != "NONE")
    {
      g4geonodename += SuperDetector();
      paramscontainer = findNode::getClass<PHParametersContainer>(parNode,g4geonodename);
      if (! paramscontainer)
	{
	  PHNodeIterator parIter(parNode);
          PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",SuperDetector()));
	  if (! DetNode)
	    {
	      DetNode = new PHCompositeNode(SuperDetector());
	      parNode->addNode(DetNode);
	    }
	  paramscontainer = new PHParametersContainer(superdetector);
	  DetNode->addNode(new PHDataNode<PHParametersContainer>(paramscontainer,g4geonodename));
	}
      paramscontainer->AddPHParameters(layer,params);
      paramnodename += superdetector;
      calibdetname = superdetector;
      isSuperDetector = 1;
    }
  else
    {
      g4geonodename += params->Name();
      parNode->addNode(new PHDataNode<PHParameters>(params,g4geonodename));
      paramnodename += params->Name();
      calibdetname = params->Name();
    }



  // ASSUMPTION: if we read from DB and/or file we don't want the stuff from
  // the node tree
  // We leave the defaults intact in case there is no entry for
  // those in the object read from the DB or file
  // Order: read first DB, then calib file if both are enabled
  if (ReadDB() || get_filetype() != PHG4DetectorSubsystem::none)
    {
      if (ReadDB())
	{
	   ReadParamsFromDB(calibdetname,isSuperDetector);
	}
      if (get_filetype() != PHG4DetectorSubsystem::none)
	{
	  ReadParamsFromFile(calibdetname, get_filetype(),isSuperDetector );
	}
    }
  else
    {
      PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode,paramnodename);
      if (nodeparams)
	{
	  params->FillFrom(nodeparams, layer);
	}
    }
  // parameters set in the macro always override whatever is read from
  // the node tree, DB or file
  UpdateParametersWithMacro();
  // save updated persistant copy on node tree
  PHCompositeNode *RunDetNode = runNode;
  if (superdetector != "NONE")
    {
      PHNodeIterator runIter(runNode);
      RunDetNode = dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",SuperDetector()));
      if (! RunDetNode)
	{
	  RunDetNode = new PHCompositeNode(SuperDetector());
	  runNode->addNode(RunDetNode);
	}
    }
  params->SaveToNodeTree(RunDetNode,paramnodename,layer);
  int iret = InitRunSubsystem(topNode);
  if (Verbosity() > 0)
    {
      PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode,paramnodename);
      cout << Name() << endl;
      nodeparams->print();
    }
  beginrunexecuted = 1;
  return iret;
}

void
PHG4DetectorSubsystem::SuperDetector(const std::string &name)
{
  superdetector = name;
  return;
}

void
PHG4DetectorSubsystem::set_double_param(const std::string &name, const double dval)
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

double
PHG4DetectorSubsystem::get_double_param(const std::string &name) const
{
  return params->get_double_param(name);
}

void
PHG4DetectorSubsystem::set_int_param(const std::string &name, const int ival)
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

int
PHG4DetectorSubsystem::get_int_param(const std::string &name) const
{
  return params->get_int_param(name);
}

void
PHG4DetectorSubsystem::set_string_param(const std::string &name, const string &sval)
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

string
PHG4DetectorSubsystem::get_string_param(const std::string &name) const
{
  return params->get_string_param(name);
}

void
PHG4DetectorSubsystem::UpdateParametersWithMacro()
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

void
PHG4DetectorSubsystem::set_default_double_param( const std::string &name, const double dval)
{
  if (default_double.find(name) == default_double.end())
    {
      default_double[name] = dval;
    }
  else
    {
      cout << "trying to overwrite default double " << name << " " 
	   << default_double[name] << " with " << dval << endl;
      exit(1);
    }
  return;
}

void
PHG4DetectorSubsystem::set_default_int_param( const std::string &name, const int ival)
{
  if (default_int.find(name) == default_int.end())
    {
      default_int[name] = ival;
    }
  else
    {
      cout << "trying to overwrite default int " << name << " " 
	   << default_int[name] << " with " << ival << endl;
      exit(1);
    }
  return;
}

void
PHG4DetectorSubsystem::set_default_string_param( const std::string &name, const string &sval)
{
  if (default_string.find(name) == default_string.end())
    {
      default_string[name] = sval;
    }
  else
    {
      cout << "trying to overwrite default string " << name << " " 
	   << default_string[name] << " with " << sval << endl;
      exit(1);
    }
  return;
}

void
PHG4DetectorSubsystem::InitializeParameters()
{
  set_default_int_param("absorberactive", 0);
  set_default_int_param("absorbertruth", 0);
  set_default_int_param("active", 0);
  set_default_int_param("blackhole", 0);

  SetDefaultParameters(); // call method from specific subsystem
  // now load those parameters to our params class
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

int
PHG4DetectorSubsystem::SaveParamsToDB()
{
  int iret = 0;
  if (paramscontainer)
    {
      iret = paramscontainer->WriteToDB();
    }
  else
    {
      iret = params->WriteToDB();
    }
  if (iret)
    {
      cout << "problem committing to DB" << endl;
    }
  return iret;
}

int
PHG4DetectorSubsystem::ReadParamsFromDB(const string &name, const int issuper)
{
  int iret = 0;
  if (issuper)
    {
      iret = params->ReadFromDB(name,layer);
    }
  else
    {
      iret = params->ReadFromDB();
    }
  if (iret)
    {
      cout << "problem reading from DB" << endl;
    }
  return iret;
}

int
PHG4DetectorSubsystem::SaveParamsToFile(const PHG4DetectorSubsystem::FILE_TYPE ftyp)
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
  int iret = 0;
  if (paramscontainer)
    {
      iret = paramscontainer->WriteToFile(extension,calibfiledir);
    }
  else
    {
      iret = params->WriteToFile(extension,calibfiledir);
    }
  if (iret)
    {
      cout << "problem saving to " << extension << " file " << endl;
    }
  return iret;
}

int
PHG4DetectorSubsystem::ReadParamsFromFile(const string &name, const PHG4DetectorSubsystem::FILE_TYPE ftyp, const int issuper)
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
  int iret = params->ReadFromFile(name, extension, layer, issuper, calibfiledir);
  if (iret)
    {
      cout << "problem reading from " << extension << " file " << endl;
    }
  return iret;
}

void
PHG4DetectorSubsystem::SetActive(const int i)
{
  iparams["active"] = i;
}

void
PHG4DetectorSubsystem::SetAbsorberActive(const int i)
{
  iparams["absorberactive"] = i;
}

void
PHG4DetectorSubsystem::BlackHole(const int i)
{
  iparams["blackhole"] = i;
}

void
PHG4DetectorSubsystem::SetAbsorberTruth(const int i)
{
  iparams["absorbertruth"] = i;
}

