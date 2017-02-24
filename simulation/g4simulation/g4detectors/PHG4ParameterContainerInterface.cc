#include "PHG4ParameterContainerInterface.h"
#include "PHG4ParametersContainer.h"
#include "PHG4Parameters.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>

#include <TSystem.h>

#include <sstream>

using namespace std;

PHG4ParameterContainerInterface::PHG4ParameterContainerInterface(const std::string &name):
  paramscontainer(new PHG4ParametersContainer(name))
{}

void
PHG4ParameterContainerInterface::set_name(const string &name)
{
  paramscontainer->set_name(name);
}

void
PHG4ParameterContainerInterface::set_default_double_param( const std::string &name, const double dval)
{
  if (default_double.find(name) == default_double.end())
    {
      default_double[name] = dval;
    }
  else
    {
      cout << "trying to overwrite default double " << name << " " 
	   << default_double[name] << " with " << dval << endl;
      gSystem->Exit(1);
    }
  return;
}

void
PHG4ParameterContainerInterface::set_default_int_param( const std::string &name, const int ival)
{
  if (default_int.find(name) == default_int.end())
    {
      default_int[name] = ival;
    }
  else
    {
      cout << "trying to overwrite default int " << name << " " 
	   << default_int[name] << " with " << ival << endl;
      gSystem->Exit(1);
    }
  return;
}

void
PHG4ParameterContainerInterface::set_default_string_param( const std::string &name, const string &sval)
{
  if (default_string.find(name) == default_string.end())
    {
      default_string[name] = sval;
    }
  else
    {
      cout << "trying to overwrite default string " << name << " " 
	   << default_string[name] << " with " << sval << endl;
      gSystem->Exit(1);
    }
  return;
}
void
PHG4ParameterContainerInterface::set_double_param(const int detid, const std::string &name, const double dval)
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
  PHG4Parameters *params = paramscontainer->GetParametersToModify(detid);
  if (! params)
    {
      ostringstream paramname;
      paramname << paramscontainer->Name() << "_" << detid;
      params = new PHG4Parameters(paramname.str());
      paramscontainer->AddPHG4Parameters(detid,params);
    }
  params->set_double_param(name,dval);
}

double
PHG4ParameterContainerInterface::get_double_param(const int detid, const std::string &name) const
{
  const PHG4Parameters *params = paramscontainer->GetParameters(detid);
  if (params)
    {
      return params->get_double_param(name);
    }
      cout << "no parameters for detid " << detid << " in " 
           << paramscontainer->Name() << " found" << endl;
  return NAN;
}

void
PHG4ParameterContainerInterface::set_int_param(const int detid, const std::string &name, const int ival)
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
  PHG4Parameters *params = paramscontainer->GetParametersToModify(detid);
  if (! params)
    {
      ostringstream paramname;
      paramname << paramscontainer->Name() << "_" << detid;
      params = new PHG4Parameters(paramname.str());
      paramscontainer->AddPHG4Parameters(detid,params);
    }
  params->set_int_param(name,ival);
}

int
PHG4ParameterContainerInterface::get_int_param(const int detid, const std::string &name) const
{
  const PHG4Parameters *params = paramscontainer->GetParameters(detid);
  if (params)
    {
      return params->get_int_param(name);
    }
      cout << "no parameters for detid " << detid << " in " 
           << paramscontainer->Name() << " found" << endl;
      return (~0x0);
}

void
PHG4ParameterContainerInterface::set_string_param(const int detid, const std::string &name, const string &sval)
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
  PHG4Parameters *params = paramscontainer->GetParametersToModify(detid);
  if (! params)
    {
      ostringstream paramname;
      paramname << paramscontainer->Name() << "_" << detid;
      params = new PHG4Parameters(paramname.str());
      paramscontainer->AddPHG4Parameters(detid,params);
    }
  params->set_string_param(name,sval);
}

string
PHG4ParameterContainerInterface::get_string_param(const int detid, const std::string &name) const
{
  const PHG4Parameters *params = paramscontainer->GetParameters(detid);
  if (params)
    {
      return params->get_string_param(name);
    }
      cout << "no parameters for detid " << detid << " in " 
           << paramscontainer->Name() << " found" << endl;
      return "";
}

void
PHG4ParameterContainerInterface::UpdateParametersWithMacro()
{
  map<int, PHG4Parameters *>::const_iterator iter;
  for (iter = macroparams.begin(); iter != macroparams.end(); ++iter)
    {
      PHG4Parameters *params = paramscontainer->GetParametersToModify(iter->first);
      std::pair< PHG4Parameters::dIter, PHG4Parameters::dIter> double_begin_end = iter->second->get_all_double_params();
      for (PHG4Parameters::dIter diter = double_begin_end.first; diter != double_begin_end.second; ++diter)
	{
	  params->set_double_param(diter->first,diter->second);
	}

      std::pair< PHG4Parameters::iIter, PHG4Parameters::iIter> int_begin_end = iter->second->get_all_int_params();
      for (PHG4Parameters::iIter iiter = int_begin_end.first; iiter != int_begin_end.second; ++iiter)
	{
	  params->set_int_param(iiter->first,iiter->second);
	}

      std::pair< PHG4Parameters::strIter, PHG4Parameters::strIter> string_begin_end = iter->second->get_all_string_params();
      for (PHG4Parameters::strIter striter = string_begin_end.first; striter != string_begin_end.second; ++striter)
	{
	  params->set_string_param(striter->first,striter->second);
	}
    }
  return;
}

void
PHG4ParameterContainerInterface::SaveToNodeTree(PHCompositeNode *runNode, const string &nodename)
{
  //  paramscontainer->SaveToNodeTree(runNode, nodename);
  return;
}

void
PHG4ParameterContainerInterface::PutOnParNode(PHCompositeNode *parNode, const string &nodename)
{
  parNode->addNode(new PHDataNode<PHG4ParametersContainer>(paramscontainer,nodename));
}

void
PHG4ParameterContainerInterface::InitializeParameters()
{
  SetDefaultParameters(); // call method from specific subsystem
}

void
PHG4ParameterContainerInterface::CreateInitialize(const int detid)
{
  PHG4Parameters *params = paramscontainer->GetParametersToModify(detid);
  if (! params)
    {
      ostringstream paramname;
      paramname << paramscontainer->Name() << "_" << detid;
      params = new PHG4Parameters(paramname.str());
      paramscontainer->AddPHG4Parameters(detid,params);
    }
  for (PHG4Parameters::dIter diter = default_double.begin(); diter != default_double.end(); ++diter)
    {
      params->set_double_param(diter->first,diter->second);
    }
  for (PHG4Parameters::iIter iiter = default_int.begin(); iiter != default_int.end(); ++iiter)
    {
      params->set_int_param(iiter->first,iiter->second);
    }
  for (PHG4Parameters::strIter striter = default_string.begin(); striter != default_string.end(); ++striter)
    {
      params->set_string_param(striter->first,striter->second);
    }
}

