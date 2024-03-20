#include "InttFeeMap.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <filesystem>  // for exists
#include <utility>     // for pair

int InttFeeMap::LoadFromFile(
    std::string const& filename)
{
  if (filename.empty())
  {
    std::cout << "int InttFeeMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout << "int InttFeeMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

int InttFeeMap::LoadFromCDB(
    std::string const& name)
{
  if (name.empty())
  {
    std::cout << "int InttFeeMap::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return 1;
  }

  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  return v_LoadFromCDBTTree(cdbttree);
}

void InttFeeMap::identify(
  std::ostream& out) const
{
  out << "InttFeeMap\n"
      << "\tBase Version\n"
      << "\tUnimplemented" << std::endl;
}

int InttFeeMap::Convert(
  InttMap::Online_s& /*unused*/,
  InttMap::Offline_s const& /*unused*/) const
{
  return 1;
}

int InttFeeMap::Convert(
  InttMap::Offline_s& /*unused*/,
  InttMap::Online_s const& /*unused*/) const
{
  return 1;
}

int InttFeeMap::Convert(
  InttMap::RawData_s& /*unused*/,
  InttMap::Online_s const& /*unused*/) const
{
  return 1;
}

int InttFeeMap::Convert(
  InttMap::Online_s& /*unused*/,
  InttMap::RawData_s const& /*unused*/) const
{
  return 1;
}

int InttFeeMap::Convert(
  InttMap::RawData_s& /*unused*/,
  InttMap::Offline_s const& /*unused*/) const
{
  return 1;
}

int InttFeeMap::Convert(
  InttMap::Offline_s& /*unused*/,
  InttMap::RawData_s const& /*unused*/) const
{
  return 1;
}

InttMap::Online_s InttFeeMap::ToOnline(
  InttMap::Offline_s const& ofl) const
{
	InttMap::Online_s onl;
	Convert(onl, ofl);
	return onl;
}

InttMap::Online_s InttFeeMap::ToOnline(
  InttMap::RawData_s const& raw) const
{
	InttMap::Online_s onl;
	Convert(onl, raw);
	return onl;
}

InttMap::RawData_s InttFeeMap::ToRawData(
  InttMap::Online_s const& onl) const
{
	InttMap::RawData_s raw;
	Convert(raw, onl);
	return raw;
}

InttMap::RawData_s InttFeeMap::ToRawData(
  InttMap::Offline_s const& ofl) const
{
	InttMap::RawData_s raw;
	Convert(raw, ofl);
	return raw;
}

InttMap::Offline_s InttFeeMap::ToOffline(
  InttMap::Online_s const& onl) const
{
	InttMap::Offline_s ofl;
	Convert(ofl, onl);
	return ofl;
}

InttMap::Offline_s InttFeeMap::ToOffline(
  InttMap::RawData_s const& raw) const
{
	InttMap::Offline_s ofl;
	Convert(ofl, raw);
	return ofl;
}

int InttFeeMap::v_LoadFromCDBTTree(
    CDBTTree& /*unused*/)
{
  return 1;
}
