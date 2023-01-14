#include "PHParameterInterface.h"
#include "PHParameters.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/phool.h>

#include <TSystem.h>

#include <iostream>
#include <utility>

PHParameterInterface::PHParameterInterface(const std::string &name)
  : m_Params(new PHParameters(name))
{
}

PHParameterInterface::~PHParameterInterface()
{
  delete m_Params;
}

void PHParameterInterface::set_paramname(const std::string &name)
{
  m_Params->set_name(name);
}

void PHParameterInterface::set_default_double_param(const std::string &name, const double dval)
{
  if (m_DefaultDoubleParMap.find(name) == m_DefaultDoubleParMap.end())
  {
    m_DefaultDoubleParMap[name] = dval;
  }
  else
  {
    std::cout << "trying to overwrite default double " << name << " "
              << m_DefaultDoubleParMap[name] << " with " << dval << std::endl;
    gSystem->Exit(1);
  }
  return;
}

void PHParameterInterface::set_default_int_param(const std::string &name, const int ival)
{
  if (m_DefaultIntParMap.find(name) == m_DefaultIntParMap.end())
  {
    m_DefaultIntParMap[name] = ival;
  }
  else
  {
    std::cout << "trying to overwrite default int " << name << " "
              << m_DefaultIntParMap[name] << " with " << ival << std::endl;
    gSystem->Exit(1);
  }
  return;
}

void PHParameterInterface::set_default_string_param(const std::string &name, const std::string &sval)
{
  if (m_DefaultStringParMap.find(name) == m_DefaultStringParMap.end())
  {
    m_DefaultStringParMap[name] = sval;
  }
  else
  {
    std::cout << "trying to overwrite default string " << name << " "
              << m_DefaultStringParMap[name] << " with " << sval << std::endl;
    gSystem->Exit(1);
  }
  return;
}
void PHParameterInterface::set_double_param(const std::string &name, const double dval)
{
  if (m_Locked)
  {
    std::cout << PHWHERE << " PHParameterInterface is locked, no modifictions allowd" << std::endl;
    gSystem->Exit(1);
  }
  if (m_DefaultDoubleParMap.find(name) == m_DefaultDoubleParMap.end())
  {
    std::cout << "double parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented double parameters are:" << std::endl;
    for (std::map<const std::string, double>::const_iterator iter = m_DefaultDoubleParMap.begin(); iter != m_DefaultDoubleParMap.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  m_DoubleParMap[name] = dval;
}

double
PHParameterInterface::get_double_param(const std::string &name) const
{
  return m_Params->get_double_param(name);
}

void PHParameterInterface::set_int_param(const std::string &name, const int ival)
{
  if (m_Locked)
  {
    std::cout << PHWHERE << " PHParameterInterface is locked, no modifictions allowd" << std::endl;
    gSystem->Exit(1);
  }
  if (m_DefaultIntParMap.find(name) == m_DefaultIntParMap.end())
  {
    std::cout << "integer parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented integer parameters are:" << std::endl;
    for (std::map<const std::string, int>::const_iterator iter = m_DefaultIntParMap.begin(); iter != m_DefaultIntParMap.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  m_IntParMap[name] = ival;
}

int PHParameterInterface::get_int_param(const std::string &name) const
{
  return m_Params->get_int_param(name);
}

void PHParameterInterface::set_string_param(const std::string &name, const std::string &sval)
{
  if (m_Locked)
  {
    std::cout << PHWHERE << " PHParameterInterface is locked, no modifictions allowd" << std::endl;
    gSystem->Exit(1);
  }
  if (m_DefaultStringParMap.find(name) == m_DefaultStringParMap.end())
  {
    std::cout << "string parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented string parameters are:" << std::endl;
    for (std::map<const std::string, std::string>::const_iterator iter = m_DefaultStringParMap.begin(); iter != m_DefaultStringParMap.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  m_StringParMap[name] = sval;
}

std::string
PHParameterInterface::get_string_param(const std::string &name) const
{
  return m_Params->get_string_param(name);
}

void PHParameterInterface::UpdateParametersWithMacro()
{
  for (std::map<const std::string, double>::const_iterator iter = m_DoubleParMap.begin(); iter != m_DoubleParMap.end(); ++iter)
  {
    m_Params->set_double_param(iter->first, iter->second);
  }
  for (std::map<const std::string, int>::const_iterator iter = m_IntParMap.begin(); iter != m_IntParMap.end(); ++iter)
  {
    m_Params->set_int_param(iter->first, iter->second);
  }
  for (std::map<const std::string, std::string>::const_iterator iter = m_StringParMap.begin(); iter != m_StringParMap.end(); ++iter)
  {
    m_Params->set_string_param(iter->first, iter->second);
  }
  return;
}

void PHParameterInterface::SaveToNodeTree(PHCompositeNode *runNode, const std::string &nodename)
{
  m_Locked = true; // no more modifications after it it on the node tree
  m_Params->SaveToNodeTree(runNode, nodename);
  return;
}

void PHParameterInterface::PutOnParNode(PHCompositeNode *parNode, const std::string &nodename)
{
  m_Locked = true; // no more modifications after it it on the node tree
  PHParameters *newparams = new PHParameters(*m_Params,m_Params->Name());
  parNode->addNode(new PHDataNode<PHParameters>(newparams, nodename));
}

void PHParameterInterface::InitializeParameters()
{
  SetDefaultParameters();  // call method from specific subsystem
  // now load those parameters to our params class
  for (std::map<const std::string, double>::const_iterator iter = m_DefaultDoubleParMap.begin(); iter != m_DefaultDoubleParMap.end(); ++iter)
  {
    m_Params->set_double_param(iter->first, iter->second);
  }
  for (std::map<const std::string, int>::const_iterator iter = m_DefaultIntParMap.begin(); iter != m_DefaultIntParMap.end(); ++iter)
  {
    m_Params->set_int_param(iter->first, iter->second);
  }
  for (std::map<const std::string, std::string>::const_iterator iter = m_DefaultStringParMap.begin(); iter != m_DefaultStringParMap.end(); ++iter)
  {
    m_Params->set_string_param(iter->first, iter->second);
  }
}

void PHParameterInterface::Print() const
{
  if (m_Params)
  {
    m_Params->Print();
  }
}
