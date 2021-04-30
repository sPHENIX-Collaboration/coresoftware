#include "PHParameterContainerInterface.h"
#include "PHParameters.h"
#include "PHParametersContainer.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>

#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

PHParameterContainerInterface::PHParameterContainerInterface(const std::string &name)
  : paramscontainer(new PHParametersContainer(name))
{
  std::string pname(name);
  pname += "_default";
  defaultparams = new PHParameters(pname);
}

PHParameterContainerInterface::~PHParameterContainerInterface()
{
  delete paramscontainer;
  delete defaultparams;
  while (macroparams.begin() != macroparams.end())
  {
    delete macroparams.begin()->second;
    macroparams.erase(macroparams.begin());
  }
  return;
}

void PHParameterContainerInterface::set_name(const std::string &name)
{
  paramscontainer->set_name(name);
}

void PHParameterContainerInterface::set_default_double_param(const std::string &name, const double dval)
{
  if (defaultparams->exist_double_param(name))
  {
    std::cout << "trying to overwrite default double " << name << " "
              << defaultparams->get_double_param(name) << " with " << dval << std::endl;
    gSystem->Exit(1);
  }
  defaultparams->set_double_param(name, dval);
  return;
}

void PHParameterContainerInterface::set_default_int_param(const std::string &name, const int ival)
{
  if (defaultparams->exist_int_param(name))
  {
    std::cout << "trying to overwrite default double " << name << " "
              << defaultparams->get_int_param(name) << " with " << ival << std::endl;
    gSystem->Exit(1);
  }
  defaultparams->set_int_param(name, ival);
  return;
}

void PHParameterContainerInterface::set_default_string_param(const std::string &name, const std::string &sval)
{
  if (defaultparams->exist_string_param(name))
  {
    std::cout << "trying to overwrite default double " << name << " "
              << defaultparams->get_string_param(name) << " with " << sval << std::endl;
    gSystem->Exit(1);
  }
  defaultparams->set_string_param(name, sval);
  return;
}

void PHParameterContainerInterface::set_double_param(const int detid, const std::string &name, const double dval)
{
  if (!defaultparams->exist_double_param(name))
  {
    std::cout << "double parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented double parameters are:" << std::endl;
    defaultparams->printdouble();
    gSystem->Exit(1);
    return;
  }
  std::map<int, PHParameters *>::iterator iter = macroparams.find(detid);
  PHParameters *params = nullptr;
  if (iter != macroparams.end())
  {
    params = iter->second;
  }
  else
  {
    std::ostringstream paramname;
    paramname << paramscontainer->Name() << "_" << detid;
    params = new PHParameters(paramname.str());
    macroparams[detid] = params;
  }
  params->set_double_param(name, dval);
}

double
PHParameterContainerInterface::get_double_param(const int detid, const std::string &name) const
{
  const PHParameters *params = paramscontainer->GetParameters(detid);
  if (params)
  {
    return params->get_double_param(name);
  }
  std::cout << "no parameters for detid " << detid << " in "
            << paramscontainer->Name() << " found" << std::endl;
  return NAN;
}

void PHParameterContainerInterface::set_int_param(const int detid, const std::string &name, const int ival)
{
  if (!defaultparams->exist_int_param(name))
  {
    std::cout << "integer parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented integer parameters are:" << std::endl;
    defaultparams->printint();
    gSystem->Exit(1);
    return;
  }
  std::map<int, PHParameters *>::iterator iter = macroparams.find(detid);
  PHParameters *params = nullptr;
  if (iter != macroparams.end())
  {
    params = iter->second;
  }
  else
  {
    std::ostringstream paramname;
    paramname << paramscontainer->Name() << "_" << detid;
    params = new PHParameters(paramname.str());
    paramscontainer->AddPHParameters(detid, params);
  }
  params->set_int_param(name, ival);
}

int PHParameterContainerInterface::get_int_param(const int detid, const std::string &name) const
{
  const PHParameters *params = paramscontainer->GetParameters(detid);
  if (params)
  {
    return params->get_int_param(name);
  }
  std::cout << "no parameters for detid " << detid << " in "
            << paramscontainer->Name() << " found" << std::endl;
  return (~0x0);
}

void PHParameterContainerInterface::set_string_param(const int detid, const std::string &name, const std::string &sval)
{
  if (!defaultparams->exist_string_param(name))
  {
    std::cout << "string parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented string parameters are:" << std::endl;
    defaultparams->printstring();
    gSystem->Exit(1);
    return;
  }

  std::map<int, PHParameters *>::iterator iter = macroparams.find(detid);
  PHParameters *params = nullptr;
  if (iter != macroparams.end())
  {
    params = iter->second;
  }
  else
  {
    std::ostringstream paramname;
    paramname << paramscontainer->Name() << "_" << detid;
    params = new PHParameters(paramname.str());
    paramscontainer->AddPHParameters(detid, params);
  }
  params->set_string_param(name, sval);
}

std::string
PHParameterContainerInterface::get_string_param(const int detid, const std::string &name) const
{
  const PHParameters *params = paramscontainer->GetParameters(detid);
  if (params)
  {
    return params->get_string_param(name);
  }
  std::cout << "no parameters for detid " << detid << " in "
            << paramscontainer->Name() << " found" << std::endl;
  return "";
}

void PHParameterContainerInterface::UpdateParametersWithMacro()
{
  std::map<int, PHParameters *>::const_iterator iter;
  for (iter = macroparams.begin(); iter != macroparams.end(); ++iter)
  {
    CreateInitialize(iter->first);

    PHParameters *params = paramscontainer->GetParametersToModify(iter->first);
    std::pair<PHParameters::dIter, PHParameters::dIter> double_begin_end = iter->second->get_all_double_params();
    for (PHParameters::dIter diter = double_begin_end.first; diter != double_begin_end.second; ++diter)
    {
      params->set_double_param(diter->first, diter->second);
    }

    std::pair<PHParameters::iIter, PHParameters::iIter> int_begin_end = iter->second->get_all_int_params();
    for (PHParameters::iIter iiter = int_begin_end.first; iiter != int_begin_end.second; ++iiter)
    {
      params->set_int_param(iiter->first, iiter->second);
    }

    std::pair<PHParameters::strIter, PHParameters::strIter> string_begin_end = iter->second->get_all_string_params();
    for (PHParameters::strIter striter = string_begin_end.first; striter != string_begin_end.second; ++striter)
    {
      params->set_string_param(striter->first, striter->second);
    }
  }
  return;
}

void PHParameterContainerInterface::SaveToNodeTree(PHCompositeNode *runNode, const std::string &nodename)
{
  paramscontainer->SaveToNodeTree(runNode, nodename);
  return;
}

void PHParameterContainerInterface::PutOnParNode(PHCompositeNode *parNode, const std::string &nodename)
{
  parNode->addNode(new PHDataNode<PHParametersContainer>(paramscontainer, nodename));
}

void PHParameterContainerInterface::InitializeParameters()
{
  SetDefaultParameters();  // call method from specific subsystem
}

void PHParameterContainerInterface::CreateInitialize(const int detid)
{
  PHParameters *params = paramscontainer->GetParametersToModify(detid);
  if (!params)
  {
    std::ostringstream paramname;
    paramname << paramscontainer->Name() << "_" << detid;
    params = new PHParameters(*defaultparams, paramname.str());
    paramscontainer->AddPHParameters(detid, params);
  }
  return;
}

int PHParameterContainerInterface::ExistDetid(const int detid) const
{
  return paramscontainer->ExistDetid(detid);
}
