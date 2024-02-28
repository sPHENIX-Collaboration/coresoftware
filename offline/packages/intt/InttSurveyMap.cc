#include "InttSurveyMap.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <filesystem>  // for exists
#include <utility>     // for pair

void InttSurveyMap::identify(
    std::ostream& out) const
{
  out << "InttSurveyMap"
      << "\n"
      << "\tBase Version"
      << "\n"
      << "\tUnimplemented" << std::endl;
}

std::size_t InttSurveyMap::size() const
{
  return 0;
}

int InttSurveyMap::LoadFromFile(
    std::string const& filename)
{
  if (filename.empty())
  {
    std::cout << "int InttSurveyMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout << "int InttSurveyMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttSurveyMap::LoadFromCDB(
    std::string const& name)
{
  if (name.empty())
  {
    std::cout << "int InttSurveyMap::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return 1;
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttSurveyMap::v_LoadFromCDBTTree(
    CDBTTree&
    /*unused*/)
{
  return 0;
}

InttSurveyMap::val_t const* InttSurveyMap::GetAbsoluteTransform(
    key_t k) const
{
  map_t::const_iterator itr;

  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  k.strip_z = InttMap::Wildcard;
  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  k.strip_phi = InttMap::Wildcard;
  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  k.ladder_z = InttMap::Wildcard;
  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  k.ladder_phi = InttMap::Wildcard;
  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  k.layer = InttMap::Wildcard;
  if (v_LookupAbsoluteTransform(k, itr))
  {
    return &(itr->second);
  }

  return (val_t const*) nullptr;
}

InttSurveyMap::val_t const* InttSurveyMap::GetRelativeTransform(
    key_t k) const
{
  map_t::const_iterator itr;

  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  k.strip_z = InttMap::Wildcard;
  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  k.strip_phi = InttMap::Wildcard;
  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  k.ladder_z = InttMap::Wildcard;
  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  k.ladder_phi = InttMap::Wildcard;
  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  k.layer = InttMap::Wildcard;
  if (v_LookupRelativeTransform(k, itr))
  {
    return &(itr->second);
  }

  return (val_t const*) nullptr;
}

int InttSurveyMap::v_LookupAbsoluteTransform(
    key_t const& /*unused*/, map_t::const_iterator&
    /*unused*/) const
{
  return 0;
}

int InttSurveyMap::v_LookupRelativeTransform(
    key_t const& /*unused*/, map_t::const_iterator&
    /*unused*/) const
{
  return 0;
}
