#include "PdbParameterMap.h"

#include <boost/functional/hash.hpp>
#include <iostream>

using namespace std;

void PdbParameterMap::print() const
{
  cout << "PdbParameterMap::print - Hash 0x" << std::hex << get_hash() << std::dec << endl;

  cout << "double parameters: " << endl;
  for (map<const string, double>::const_iterator iter = dparams.begin(); iter != dparams.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
  }
  cout << "integer parameters: " << endl;
  for (map<const string, int>::const_iterator iter = iparams.begin(); iter != iparams.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
  }
  cout << "string parameters: " << endl;
  for (map<const string, string>::const_iterator iter = cparams.begin(); iter != cparams.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
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

  for (dMap::const_iterator iter = dparams.begin();
       iter != dparams.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (iMap::const_iterator iter = iparams.begin();
       iter != iparams.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (strMap::const_iterator iter = cparams.begin();
       iter != cparams.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  return seed;
}
