#include "PdbParameterMap.h"

#include <boost/functional/hash.hpp>

#include <iostream>

void PdbParameterMap::print() const
{
  std::cout << "PdbParameterMap::print - Hash 0x" << std::hex << get_hash() << std::dec << std::endl;

  std::cout << "double parameters: " << std::endl;
  for (const auto &dparam : dparams)
  {
    std::cout << dparam.first << ": " << dparam.second << std::endl;
  }
  std::cout << "integer parameters: " << std::endl;
  for (const auto &iparam : iparams)
  {
    std::cout << iparam.first << ": " << iparam.second << std::endl;
  }
  std::cout << "string parameters: " << std::endl;
  for (const auto &cparam : cparams)
  {
    std::cout << cparam.first << ": " << cparam.second << std::endl;
  }
}

void PdbParameterMap::Reset()
{
  dparams.clear();
  iparams.clear();
  cparams.clear();
  return;
}

void PdbParameterMap::set_int_param(const std::string &name, const int ival)
{
  iparams[name] = ival;
}

void PdbParameterMap::set_double_param(const std::string &name, const double dval)
{
  dparams[name] = dval;
}

void PdbParameterMap::set_string_param(const std::string &name, const std::string &str)
{
  cparams[name] = str;
}

size_t
PdbParameterMap::get_hash() const
{
  size_t seed = 0;

  for (const auto &dparam : dparams)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, dparam.first);
    boost::hash_combine(seed, dparam.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  for (const auto &iparam : iparams)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iparam.first);
    boost::hash_combine(seed, iparam.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  for (const auto &cparam : cparams)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, cparam.first);
    boost::hash_combine(seed, cparam.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  return seed;
}
