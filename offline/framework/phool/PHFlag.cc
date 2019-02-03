#include "PHFlag.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

const string
PHFlag::get_CharFlag(const string &name) const
{
  map<string, string>::const_iterator iter = charflag.find(name);
  if (iter != charflag.end())
  {
    return iter->second;
  }
  cout << "PHFlag::getString: ERROR Unknown character Flag " << name
       << ", The following are implemented: " << endl;
  Print();
  return nullptr;
}

const string
PHFlag::get_CharFlag(const string &name, const string &defaultval)
{
  map<string, string>::const_iterator iter = charflag.find(name);
  if (iter != charflag.end())
  {
    return iter->second;
  }
  else
  {
    set_CharFlag(name, defaultval);
    return get_CharFlag(name);
  }
}

void PHFlag::set_CharFlag(const string &name, const string &charstr)
{
  charflag[name] = charstr;
  return;
}

double PHFlag::get_DoubleFlag(const string &name) const
{
  map<string, double>::const_iterator iter = doubleflag.find(name);
  if (iter != doubleflag.end())
  {
    return iter->second;
  }
  cout << "PHFlag::getFlag: ERROR Unknown Double Flag " << name
       << ", The following are implemented: " << endl;
  Print();
  return 0.0;
}

double PHFlag::get_DoubleFlag(const string &name, const double defaultval)
{
  map<string, double>::const_iterator iter = doubleflag.find(name);
  if (iter != doubleflag.end())
  {
    return iter->second;
  }
  else
  {
    set_DoubleFlag(name, defaultval);
    return get_DoubleFlag(name);
  }
}

void PHFlag::set_DoubleFlag(const string &name, const double iflag)
{
  doubleflag[name] = iflag;
  return;
}

float PHFlag::get_FloatFlag(const string &name) const
{
  map<string, float>::const_iterator iter = floatflag.find(name);
  if (iter != floatflag.end())
  {
    return iter->second;
  }
  cout << "PHFlag::getFlag: ERROR Unknown Float Flag " << name
       << ", The following are implemented: " << endl;
  Print();
  return 0.0;
}

float PHFlag::get_FloatFlag(const string &name, const float defaultval)
{
  map<string, float>::const_iterator iter = floatflag.find(name);
  if (iter != floatflag.end())
  {
    return iter->second;
  }
  else
  {
    set_FloatFlag(name, defaultval);
    return get_FloatFlag(name);
  }
}

void PHFlag::set_FloatFlag(const string &name, const float iflag)
{
  floatflag[name] = iflag;
  return;
}

int PHFlag::get_IntFlag(const string &name) const
{
  map<string, int>::const_iterator iter = intflag.find(name);
  if (iter != intflag.end())
  {
    return iter->second;
  }
  cout << "PHFlag::getFlag: ERROR Unknown Int Flag " << name
       << ", The following are implemented: " << endl;
  Print();
  return 0;
}

int PHFlag::get_IntFlag(const string &name, int defaultval)
{
  map<string, int>::const_iterator iter = intflag.find(name);
  if (iter != intflag.end())
  {
    return iter->second;
  }
  else
  {
    set_IntFlag(name, defaultval);
    return get_IntFlag(name);
  }
}

void PHFlag::set_IntFlag(const string &name, const int iflag)
{
  intflag[name] = iflag;
  return;
}

void PHFlag::Print() const
{
  PrintIntFlags();
  PrintFloatFlags();
  PrintDoubleFlags();
  PrintCharFlags();
  return;
}

void PHFlag::PrintIntFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  cout << endl
       << "Integer Flags:" << endl;
  map<string, int>::const_iterator intiter;
  for (intiter = intflag.begin(); intiter != intflag.end(); ++intiter)
  {
    cout << intiter->first << " is " << intiter->second << endl;
  }
  return;
}

void PHFlag::PrintDoubleFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  cout << endl
       << "Double Flags:" << endl;
  map<string, double>::const_iterator doubleiter;
  for (doubleiter = doubleflag.begin(); doubleiter != doubleflag.end(); ++doubleiter)
  {
    cout << doubleiter->first << " is " << doubleiter->second << endl;
  }
  return;
}

void PHFlag::PrintFloatFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  cout << endl
       << "Float Flags:" << endl;
  map<string, float>::const_iterator floatiter;
  for (floatiter = floatflag.begin(); floatiter != floatflag.end(); ++floatiter)
  {
    cout << floatiter->first << " is " << floatiter->second << endl;
  }
  return;
}

void PHFlag::PrintCharFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  cout << endl
       << "char* Flags:" << endl;
  map<string, string>::const_iterator chariter;
  for (chariter = charflag.begin(); chariter != charflag.end(); ++chariter)
  {
    cout << chariter->first << " is " << chariter->second << endl;
  }
  return;
}

int PHFlag::FlagExist(const string &name) const
{
  map<string, int>::const_iterator iter = intflag.find(name);
  if (iter != intflag.end())
  {
    return 1;
  }
  map<string, float>::const_iterator fiter = floatflag.find(name);
  if (fiter != floatflag.end())
  {
    return 1;
  }
  map<string, double>::const_iterator diter = doubleflag.find(name);
  if (diter != doubleflag.end())
  {
    return 1;
  }
  map<string, string>::const_iterator citer = charflag.find(name);
  if (citer != charflag.end())
  {
    return 1;
  }
  return 0;
}

void PHFlag::ReadFromFile(const string &name)
{
  string label;
  float fvalue;
  int fvaluecount = 0;
  double dvalue;
  int dvaluecount = 0;
  int ivalue;
  int ivaluecount = 0;
  string cvalue;
  int cvaluecount = 0;
  string junk;
  int junkcount = 0;

  ifstream infile(name.c_str());
  while (infile >> label)
  {
    cout << "Label" << label;
    if (label.substr(0, 1) == "C")
    {
      infile >> cvalue;
      cvaluecount++;
      set_CharFlag(label.substr(1, label.size() - 1), cvalue);
      cout << " C read " << cvalue << endl;
    }
    else if (label.substr(0, 1) == "F")
    {
      infile >> fvalue;
      fvaluecount++;
      set_FloatFlag(label.substr(1, label.size() - 1), fvalue);
      cout << " F read " << fvalue << endl;
    }
    else if (label.substr(0, 1) == "D")
    {
      infile >> dvalue;
      dvaluecount++;
      set_DoubleFlag(label.substr(1, label.size() - 1), dvalue);
      cout << " D read " << dvalue << endl;
    }
    else if (label.substr(0, 1) == "I")
    {
      infile >> ivalue;
      ivaluecount++;
      set_IntFlag(label.substr(1, label.size() - 1), ivalue);
      cout << " I read " << ivalue << endl;
    }
    else
    {
      infile >> junk;
      junkcount++;
      cout << " Junk read " << junk << endl;
    }
  }

  cout << "Read CharFlags(" << cvaluecount
       << ") FloatFlags(" << fvaluecount
       << ") DoubleFlags(" << dvaluecount
       << ") IntFlags(" << ivaluecount
       << ") JunkEntries(" << junkcount
       << ") from file " << name << endl;

  infile.close();
}

void PHFlag::WriteToFile(const string &name)
{
  ofstream outFile(name.c_str());
  // loop over the map and write out the content
  map<string, int>::const_iterator intiter;
  for (intiter = intflag.begin(); intiter != intflag.end(); ++intiter)
  {
    outFile << "I" << intiter->first << "\t" << intiter->second << endl;
  }

  map<string, float>::const_iterator floatiter;
  for (floatiter = floatflag.begin(); floatiter != floatflag.end(); ++floatiter)
  {
    outFile << "F" << floatiter->first << "\t" << floatiter->second << endl;
  }

  int oldprecision = outFile.precision(15);
  map<string, double>::const_iterator doubleiter;
  for (doubleiter = doubleflag.begin(); doubleiter != doubleflag.end(); ++doubleiter)
  {
    outFile << "D" << doubleiter->first << "\t" << doubleiter->second << endl;
  }
  outFile.precision(oldprecision);

  map<string, string>::const_iterator chariter;
  for (chariter = charflag.begin(); chariter != charflag.end(); ++chariter)
  {
    outFile << "C" << chariter->first << "\t" << chariter->second << endl;
  }

  outFile.close();
}
