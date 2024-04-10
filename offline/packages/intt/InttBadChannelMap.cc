#include "InttBadChannelMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <filesystem>  // for exists

int InttBadChannelMap::LoadFromFile(
    std::string const& filename)
{
  if (filename.empty())
  {
    std::cout << "int InttBadChannelMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout << "int InttBadChannelMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttBadChannelMap::LoadFromCDB(
    std::string const& name)
{
  if (name.empty())
  {
    std::cout << "int InttBadChannelMap::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return 1;
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

void InttBadChannelMap::identify(
    std::ostream& out) const
{
  out << "InttBadChannelMap\n"
      << "\tBase Version\n"
      << "\tUnimplemented" << std::endl;
}

std::size_t InttBadChannelMap::size() const
{
  return 0;
}

bool InttBadChannelMap::IsBad(
    InttMap::Online_s const&
    /*unused*/) const
{
  std::cout << "InttBadChannelMap::IsBad\n"
           << "\tUnplemented overload (InttMap::Online_s)" << std::endl;
  return false;
}

bool InttBadChannelMap::IsBad(
    InttMap::Offline_s const&
    /*unused*/) const
{
  std::cout << "InttBadChannelMap::IsBad\n"
           << "\tUnplemented overload (InttMap::Offline_s)" << std::endl;
  return false;
}

bool InttBadChannelMap::IsBad(
    InttMap::RawData_s const&
    /*unused*/) const
{
  std::cout << "InttBadChannelMap::IsBad\n"
           << "\tUnplemented overload (InttMap::RawData_s)" << std::endl;
  return false;
}

int InttBadChannelMap::v_LoadFromCDBTTree(
    CDBTTree&
    /*unused*/)
{
  return 0;
}
