#include "InttLoadable.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <iostream>
#include <filesystem>

int InttLoadable::Load(
    std::string name)
{
  m_loaded = 0;

  if (name.empty())
  {
    name = DefaultName();
  }

  if(name.find(".root") == std::string::npos)
  {
    name = CDBInterface::instance()->getUrl(name);
    if (name.empty())
    {
      std::cerr << __PRETTY_FUNCTION__ << std::endl;
      std::cerr << "\tPayload type \"" << DefaultName() << "\" not found in CDB" << std::endl;
      return 1;
    }
  }
  else if (!std::filesystem::exists(name))
  {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tFile \"" << name << "\" does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(name);
  cdbttree.LoadCalibrations();

  m_loaded = (LoadFromCDBTTree(cdbttree) == 0);
  return m_loaded ? 0 : 1;
}

bool InttLoadable::IsLoaded() const
{
  return m_loaded;
}

int InttLoadable::Verbosity() const
{
  return m_verbosity;
}

int InttLoadable::Verbosity(
  int const& v
)
{
  return m_verbosity = v;
}

int InttLoadable::LoadFromCDBTTree(
  CDBTTree& /* unused */ )
{
  return 1;
}

