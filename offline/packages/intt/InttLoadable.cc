#include "InttLoadable.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <iostream>
#include <filesystem>

int InttLoadable::LoadFromFile(
    std::string const& filename)
{
  if (filename.empty())
  {
    return LoadFromFile();
  }

  if (!std::filesystem::exists(filename))
  {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return LoadFromCDBTTree(cdbttree);
}

int InttLoadable::LoadFromCDB(
    std::string const& name)
{
  if (name.empty())
  {
    return LoadFromCDB();
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return LoadFromCDBTTree(cdbttree);
}

int InttLoadable::LoadFromFile()
{
  std::cerr << __PRETTY_FUNCTION__ << "\n"
            << "\tDynamic function call was to an instance of base class\n"
            << "\tThis virtual function should be overridden by children\n"
            << "\tto call InttLoadable::LoadFromFile(std::string const&)\n"
            << "\twith a befitting default argument, creating the effect\n"
            << "\tof an having an overloaded method with a default value" << std::endl;
  return 1;
}

int InttLoadable::LoadFromCDB()
{
  std::cerr << __PRETTY_FUNCTION__ << "\n"
            << "\tDynamic function call was to an instance of base class\n"
            << "\tThis virtual function should be overridden by children\n"
            << "\tto call InttLoadable::LoadFromFile(std::string const&)\n"
            << "\twith a befitting default argument, creating the effect\n"
            << "\tof an having an overloaded method with a default value" << std::endl;
  return 1;
}

int InttLoadable::LoadFromCDBTTree(
    CDBTTree& /*unused*/)
{
  return 0;
}
