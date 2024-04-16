#include "InttLoadable.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <iostream>
#include <filesystem>

int InttLoadable::LoadFromFile(
    std::string filename)
{
  m_loaded = 0;
  if (filename.empty())
  {
    filename = DefaultFileName();
  }

  if (!std::filesystem::exists(filename))
  {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tFile \"" << filename << "\" does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  m_loaded = (LoadFromCDBTTree(cdbttree) == 0);
  return m_loaded ? 0 : 1;
}

int InttLoadable::LoadFromCDB(
    std::string name)
{
  m_loaded = 0;

  if (name.empty())
  {
    name = DefaultCDBName();
  }
  name = CDBInterface::instance()->getUrl(name);
  if (name.empty())
  {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tPayload type \"" << DefaultCDBName() << "\" not found in CDB" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(name);
  cdbttree.LoadCalibrations();

  m_loaded = (LoadFromCDBTTree(cdbttree) == 0);
  return m_loaded ? 0 : 1;
}

int InttLoadable::Loaded() const
{
  return m_loaded;
}

