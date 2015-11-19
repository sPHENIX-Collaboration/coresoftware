#include "PdbParameterMap.h"

#include <iostream>

using namespace std;

void
PdbParameterMap::print() const
{
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

void
PdbParameterMap::set_int_param(const std::string &name, const int ival)
{
  iparams[name] = ival;
}

void
PdbParameterMap::set_double_param(const std::string &name, const double dval)
{
  dparams[name] = dval;
}

void
PdbParameterMap::set_string_param(const std::string &name, const string &str)
{
  cparams[name] = str;
}

