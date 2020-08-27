#include "PHParameterInterface.h"
#include "PHParameters.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>

#include <TSystem.h>

#include <iostream>
#include <utility>

using namespace std;

PHParameterInterface::PHParameterInterface(const std::string &name)
  : m_Params(new PHParameters(name))
{
}

void PHParameterInterface::set_paramname(const string &name)
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
    cout << "trying to overwrite default double " << name << " "
         << m_DefaultDoubleParMap[name] << " with " << dval << endl;
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
    cout << "trying to overwrite default int " << name << " "
         << m_DefaultIntParMap[name] << " with " << ival << endl;
    gSystem->Exit(1);
  }
  return;
}

void PHParameterInterface::set_default_string_param(const std::string &name, const string &sval)
{
  if (m_DefaultStringParMap.find(name) == m_DefaultStringParMap.end())
  {
    m_DefaultStringParMap[name] = sval;
  }
  else
  {
    cout << "trying to overwrite default string " << name << " "
         << m_DefaultStringParMap[name] << " with " << sval << endl;
    gSystem->Exit(1);
  }
  return;
}
void PHParameterInterface::set_double_param(const std::string &name, const double dval)
{
  if (m_DefaultDoubleParMap.find(name) == m_DefaultDoubleParMap.end())
  {
    cout << "double parameter " << name << " not implemented" << endl;
    cout << "implemented double parameters are:" << endl;
    for (map<const string, double>::const_iterator iter = m_DefaultDoubleParMap.begin(); iter != m_DefaultDoubleParMap.end(); ++iter)
    {
      cout << iter->first << endl;
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
  if (m_DefaultIntParMap.find(name) == m_DefaultIntParMap.end())
  {
    cout << "integer parameter " << name << " not implemented" << endl;
    cout << "implemented integer parameters are:" << endl;
    for (map<const string, int>::const_iterator iter = m_DefaultIntParMap.begin(); iter != m_DefaultIntParMap.end(); ++iter)
    {
      cout << iter->first << endl;
    }
    return;
  }
  m_IntParMap[name] = ival;
}

int PHParameterInterface::get_int_param(const std::string &name) const
{
  return m_Params->get_int_param(name);
}

void PHParameterInterface::set_string_param(const std::string &name, const string &sval)
{
  if (m_DefaultStringParMap.find(name) == m_DefaultStringParMap.end())
  {
    cout << "string parameter " << name << " not implemented" << endl;
    cout << "implemented string parameters are:" << endl;
    for (map<const string, string>::const_iterator iter = m_DefaultStringParMap.begin(); iter != m_DefaultStringParMap.end(); ++iter)
    {
      cout << iter->first << endl;
    }
    return;
  }
  m_StringParMap[name] = sval;
}

string
PHParameterInterface::get_string_param(const std::string &name) const
{
  return m_Params->get_string_param(name);
}

void PHParameterInterface::UpdateParametersWithMacro()
{
  for (map<const string, double>::const_iterator iter = m_DoubleParMap.begin(); iter != m_DoubleParMap.end(); ++iter)
  {
    m_Params->set_double_param(iter->first, iter->second);
  }
  for (map<const string, int>::const_iterator iter = m_IntParMap.begin(); iter != m_IntParMap.end(); ++iter)
  {
    m_Params->set_int_param(iter->first, iter->second);
  }
  for (map<const string, string>::const_iterator iter = m_StringParMap.begin(); iter != m_StringParMap.end(); ++iter)
  {
    m_Params->set_string_param(iter->first, iter->second);
  }
  return;
}

void PHParameterInterface::SaveToNodeTree(PHCompositeNode *runNode, const string &nodename)
{
  m_Params->SaveToNodeTree(runNode, nodename);
  return;
}

void PHParameterInterface::PutOnParNode(PHCompositeNode *parNode, const string &nodename)
{
  parNode->addNode(new PHDataNode<PHParameters>(m_Params, nodename));
}

void PHParameterInterface::InitializeParameters()
{
  SetDefaultParameters();  // call method from specific subsystem
  // now load those parameters to our params class
  for (map<const string, double>::const_iterator iter = m_DefaultDoubleParMap.begin(); iter != m_DefaultDoubleParMap.end(); ++iter)
  {
    m_Params->set_double_param(iter->first, iter->second);
  }
  for (map<const string, int>::const_iterator iter = m_DefaultIntParMap.begin(); iter != m_DefaultIntParMap.end(); ++iter)
  {
    m_Params->set_int_param(iter->first, iter->second);
  }
  for (map<const string, string>::const_iterator iter = m_DefaultStringParMap.begin(); iter != m_DefaultStringParMap.end(); ++iter)
  {
    m_Params->set_string_param(iter->first, iter->second);
  }
}
