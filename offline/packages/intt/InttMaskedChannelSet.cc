#include "InttMaskedChannelSet.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <filesystem>  // for exists

void InttMaskedChannelSet::identify(
    std::ostream& out) const
{
  out << "InttMaskedChannelSet"
      << "\n"
      << "\tBase Version"
      << "\n"
      << "\tUnimplemented" << std::endl;
}

std::size_t InttMaskedChannelSet::size() const
{
  return 0;
}

int InttMaskedChannelSet::LoadFromFile(
    std::string const& filename)
{
  if (filename.empty())
  {
    std::cout << "int InttMaskedChannelSet::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout << "int InttMaskedChannelSet::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttMaskedChannelSet::LoadFromCDB(
    std::string const& name)
{
  if (name.empty())
  {
    std::cout << "int InttMaskedChannelSet::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return 1;
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

bool InttMaskedChannelSet::IsDeadChannel(
    int const& layer,
    int const& ladder_phi,
    int const& ladder_z,
    int const& strip_z,
    int const& strip_phi) const
{
  return IsDeadChannel((InttMap::Offline_s){
      .layer = layer,
      .ladder_phi = ladder_phi,
      .ladder_z = ladder_z,
      .strip_phi = strip_phi,
      .strip_z = strip_z,
  });
}

bool InttMaskedChannelSet::IsDeadChannel(
    InttMap::Offline_s const&
    /*unused*/) const
{
  return false;
}

int InttMaskedChannelSet::v_LoadFromCDBTTree(
    CDBTTree&
    /*unused*/)
{
  return 0;
}
