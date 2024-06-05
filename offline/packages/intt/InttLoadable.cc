#include "InttLoadable.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <iostream>
#include <filesystem>  // for exists

int InttLoadable::Load(std::string const& name)
{
  m_is_loaded = false;

  if (name.empty())
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\tArgument 'name' is empty string" << std::endl;
	return 1;
  }

  std::string filename = name.find(".root") != std::string::npos ? name : CDBInterface::instance()->getUrl(name);
  if (filename.empty())
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\tTag '" << name << "' not found in CDB" << std::endl;
	return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  m_is_loaded = (LoadFromCdbTTree(cdbttree) == 0); // 0 on success, 1 on failure
  return m_is_loaded ? 0 : 1;                      // 0 on success, 1 on failure
}

int InttLoadable::LoadFromCdbTTree(CDBTTree&)
{
  std::cerr << __PRETTY_FUNCTION__ << "\n"
            << "\tCall to unimplemented base" << std::endl;
  return 1;
}
