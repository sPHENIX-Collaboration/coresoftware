#include "PHG4ParameterContainerInterface.h"
#include "PHG4ParametersContainer.h"
#include "PHG4Parameters.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>

#include <TSystem.h>

#include <sstream>

using namespace std;

PHG4ParameterContainerInterface::PHG4ParameterContainerInterface(const std::string &name):
  paramscontainer(new PHG4ParametersContainer(name)),
  defaultparams(nullptr)
{
  string pname(name);
  pname += "_default";
  defaultparams = new PHG4Parameters(pname);
}

PHG4ParameterContainerInterface::~PHG4ParameterContainerInterface()
{
  delete paramscontainer;
  delete defaultparams;
  while(macroparams.begin() != macroparams.end())
    {
      delete macroparams.begin()->second;
      macroparams.erase(macroparams.begin());
    }
  return;
}


void
PHG4ParameterContainerInterface::set_name(const string &name)
{
  paramscontainer->set_name(name);
}

void
PHG4ParameterContainerInterface::set_default_double_param( const std::string &name, const double dval)
{
  if (defaultparams->exist_double_param(name))
    {
      cout << "trying to overwrite default double " << name << " " 
	   << defaultparams->get_double_param(name) << " with " << dval << endl;
      gSystem->Exit(1);
    }
  defaultparams->set_double_param(name,dval);
  return;
}

void
PHG4ParameterContainerInterface::set_default_int_param( const std::string &name, const int ival)
{
  if (defaultparams->exist_int_param(name))
    {
      cout << "trying to overwrite default double " << name << " " 
	   << defaultparams->get_int_param(name) << " with " << ival << endl;
      gSystem->Exit(1);
    }
  defaultparams->set_int_param(name,ival);
  return;
}

void
PHG4ParameterContainerInterface::set_default_string_param( const std::string &name, const string &sval)
{
  if (defaultparams->exist_string_param(name))
    {
      cout << "trying to overwrite default double " << name << " " 
	   << defaultparams->get_string_param(name) << " with " << sval << endl;
      gSystem->Exit(1);
    }
  defaultparams->set_string_param(name,sval);
  return;
}

void
PHG4ParameterContainerInterface::set_double_param(const int detid, const std::string &name, const double dval)
{
  if (! defaultparams->exist_double_param(name))
    {
      cout << "double parameter " << name << " not implemented" << endl;
      cout << "implemented double parameters are:" << endl;
      defaultparams->printdouble();
      gSystem->Exit(1);
      return;
    }
  map<int, PHG4Parameters *>::iterator iter = macroparams.find(detid);
  PHG4Parameters *params = nullptr;
  if (iter != macroparams.end())
    {
      params = iter->second;
    }
  else
    {
      ostringstream paramname;
      paramname << paramscontainer->Name() << "_" << detid;
      params = new PHG4Parameters(paramname.str());
      macroparams[detid] = params;
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
  if (! defaultparams->exist_int_param(name))
    {
      cout << "integer parameter " << name << " not implemented" << endl;
      cout << "implemented integer parameters are:" << endl;
      defaultparams->printint();
      gSystem->Exit(1);
      return;
    }
  map<int, PHG4Parameters *>::iterator iter = macroparams.find(detid);
  PHG4Parameters *params = nullptr;
  if (iter != macroparams.end())
    {
      params = iter->second;
    }
  else
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
  if (! defaultparams->exist_string_param(name))
    {
      cout << "string parameter " << name << " not implemented" << endl;
      cout << "implemented string parameters are:" << endl;
      defaultparams->printstring();
      gSystem->Exit(1);
      return;
    }

  map<int, PHG4Parameters *>::iterator iter = macroparams.find(detid);
  PHG4Parameters *params = nullptr;
  if (iter != macroparams.end())
    {
      params = iter->second;
    }
  else
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
      CreateInitialize(iter->first);

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
  paramscontainer->SaveToNodeTree(runNode, nodename);
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
      params = new PHG4Parameters(*defaultparams,paramname.str());
      paramscontainer->AddPHG4Parameters(detid,params);
    }
  return;
}

int
PHG4ParameterContainerInterface::ExistDetid(const int detid) const
{
  return paramscontainer->ExistDetid(detid);
}
