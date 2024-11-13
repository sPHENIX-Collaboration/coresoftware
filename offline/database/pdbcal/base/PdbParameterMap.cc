#include "PdbParameterMap.h"

#include <boost/functional/hash.hpp>
#include <iostream>

using namespace std;

void PdbParameterMap::print() const
{
  cout << "PdbParameterMap::print - Hash 0x" << std::hex << get_hash() << std::dec << endl;

  cout << "double parameters: " << endl;
  for (const auto &dparam : dparams)
  {
    cout << dparam.first << ": " << dparam.second << endl;
  }
  cout << "integer parameters: " << endl;
  for (const auto &iparam : iparams)
  {
    cout << iparam.first << ": " << iparam.second << endl;
  }
  cout << "string parameters: " << endl;
  for (const auto &cparam : cparams)
  {
    cout << cparam.first << ": " << cparam.second << endl;
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

void PdbParameterMap::set_string_param(const std::string &name, const string &str)
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
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (const auto &iparam : iparams)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iparam.first);
    boost::hash_combine(seed, iparam.second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (const auto &cparam : cparams)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, cparam.first);
    boost::hash_combine(seed, cparam.second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  return seed;
}
