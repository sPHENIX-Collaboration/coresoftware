#include "PHFlag.h"

// boost stacktrace header causes a shadow warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <fstream>
#include <iostream>
#include <map>
#include <utility>  // for pair

const std::string
PHFlag::get_CharFlag(const std::string &name) const
{
  std::map<std::string, std::string>::const_iterator iter = charflag.find(name);
  if (iter != charflag.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getString: ERROR Unknown character Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return nullptr;
}

const std::string
PHFlag::get_CharFlag(const std::string &name, const std::string &defaultval)
{
  std::map<std::string, std::string>::const_iterator iter = charflag.find(name);
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

void PHFlag::set_CharFlag(const std::string &name, const std::string &charstr)
{
  charflag[name] = charstr;
  return;
}

double PHFlag::get_DoubleFlag(const std::string &name) const
{
  std::map<std::string, double>::const_iterator iter = doubleflag.find(name);
  if (iter != doubleflag.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getFlag: ERROR Unknown Double Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return 0.0;
}

double PHFlag::get_DoubleFlag(const std::string &name, const double defaultval)
{
  std::map<std::string, double>::const_iterator iter = doubleflag.find(name);
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

void PHFlag::set_DoubleFlag(const std::string &name, const double iflag)
{
  doubleflag[name] = iflag;
  return;
}

float PHFlag::get_FloatFlag(const std::string &name) const
{
  std::map<std::string, float>::const_iterator iter = floatflag.find(name);
  if (iter != floatflag.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getFlag: ERROR Unknown Float Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return 0.0;
}

float PHFlag::get_FloatFlag(const std::string &name, const float defaultval)
{
  std::map<std::string, float>::const_iterator iter = floatflag.find(name);
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

void PHFlag::set_FloatFlag(const std::string &name, const float iflag)
{
  floatflag[name] = iflag;
  return;
}

int PHFlag::get_IntFlag(const std::string &name) const
{
  std::map<std::string, int>::const_iterator iter = intflag.find(name);
  if (iter != intflag.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getFlag: ERROR Unknown Int Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return 0;
}

int PHFlag::get_IntFlag(const std::string &name, int defaultval)
{
  std::map<std::string, int>::const_iterator iter = intflag.find(name);
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

void PHFlag::set_IntFlag(const std::string &name, const int iflag)
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
  std::cout << std::endl
       << "Integer Flags:" << std::endl;
  std::map<std::string, int>::const_iterator intiter;
  for (intiter = intflag.begin(); intiter != intflag.end(); ++intiter)
  {
    std::cout << intiter->first << " is " << intiter->second << std::endl;
  }
  return;
}

void PHFlag::PrintDoubleFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
       << "Double Flags:" << std::endl;
  std::map<std::string, double>::const_iterator doubleiter;
  for (doubleiter = doubleflag.begin(); doubleiter != doubleflag.end(); ++doubleiter)
  {
    std::cout << doubleiter->first << " is " << doubleiter->second << std::endl;
  }
  return;
}

void PHFlag::PrintFloatFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
       << "Float Flags:" << std::endl;
  std::map<std::string, float>::const_iterator floatiter;
  for (floatiter = floatflag.begin(); floatiter != floatflag.end(); ++floatiter)
  {
    std::cout << floatiter->first << " is " << floatiter->second << std::endl;
  }
  return;
}

void PHFlag::PrintCharFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
       << "char* Flags:" << std::endl;
  std::map<std::string, std::string>::const_iterator chariter;
  for (chariter = charflag.begin(); chariter != charflag.end(); ++chariter)
  {
    std::cout << chariter->first << " is " << chariter->second << std::endl;
  }
  return;
}

int PHFlag::FlagExist(const std::string &name) const
{
  std::map<std::string, int>::const_iterator iter = intflag.find(name);
  if (iter != intflag.end())
  {
    return 1;
  }
  std::map<std::string, float>::const_iterator fiter = floatflag.find(name);
  if (fiter != floatflag.end())
  {
    return 1;
  }
  std::map<std::string, double>::const_iterator diter = doubleflag.find(name);
  if (diter != doubleflag.end())
  {
    return 1;
  }
  std::map<std::string, std::string>::const_iterator citer = charflag.find(name);
  if (citer != charflag.end())
  {
    return 1;
  }
  return 0;
}

void PHFlag::ReadFromFile(const std::string &name)
{
  std::string label;
  float fvalue;
  int fvaluecount = 0;
  double dvalue;
  int dvaluecount = 0;
  int ivalue;
  int ivaluecount = 0;
  std::string cvalue;
  int cvaluecount = 0;
  std::string junk;
  int junkcount = 0;

  std::ifstream infile(name);
  while (infile >> label)
  {
    std::cout << "Label" << label;
    if (label.substr(0, 1) == "C")
    {
      infile >> cvalue;
      cvaluecount++;
      set_CharFlag(label.substr(1, label.size() - 1), cvalue);
      std::cout << " C read " << cvalue << std::endl;
    }
    else if (label.substr(0, 1) == "F")
    {
      infile >> fvalue;
      fvaluecount++;
      set_FloatFlag(label.substr(1, label.size() - 1), fvalue);
      std::cout << " F read " << fvalue << std::endl;
    }
    else if (label.substr(0, 1) == "D")
    {
      infile >> dvalue;
      dvaluecount++;
      set_DoubleFlag(label.substr(1, label.size() - 1), dvalue);
      std::cout << " D read " << dvalue << std::endl;
    }
    else if (label.substr(0, 1) == "I")
    {
      infile >> ivalue;
      ivaluecount++;
      set_IntFlag(label.substr(1, label.size() - 1), ivalue);
      std::cout << " I read " << ivalue << std::endl;
    }
    else
    {
      infile >> junk;
      junkcount++;
      std::cout << " Junk read " << junk << std::endl;
    }
  }

  std::cout << "Read CharFlags(" << cvaluecount
       << ") FloatFlags(" << fvaluecount
       << ") DoubleFlags(" << dvaluecount
       << ") IntFlags(" << ivaluecount
       << ") JunkEntries(" << junkcount
       << ") from file " << name << std::endl;

  infile.close();
}

void PHFlag::WriteToFile(const std::string &name)
{
  std::ofstream outFile(name);
  // loop over the map and write out the content
  std::map<std::string, int>::const_iterator intiter;
  for (intiter = intflag.begin(); intiter != intflag.end(); ++intiter)
  {
    outFile << "I" << intiter->first << "\t" << intiter->second << std::endl;
  }

  std::map<std::string, float>::const_iterator floatiter;
  for (floatiter = floatflag.begin(); floatiter != floatflag.end(); ++floatiter)
  {
    outFile << "F" << floatiter->first << "\t" << floatiter->second << std::endl;
  }

  int oldprecision = outFile.precision(15);
  std::map<std::string, double>::const_iterator doubleiter;
  for (doubleiter = doubleflag.begin(); doubleiter != doubleflag.end(); ++doubleiter)
  {
    outFile << "D" << doubleiter->first << "\t" << doubleiter->second << std::endl;
  }
  outFile.precision(oldprecision);

  std::map<std::string, std::string>::const_iterator chariter;
  for (chariter = charflag.begin(); chariter != charflag.end(); ++chariter)
  {
    outFile << "C" << chariter->first << "\t" << chariter->second << std::endl;
  }

  outFile.close();
}

void PHFlag::PrintStackTrace() const
{
  std::cout << "Called by #3 or #4 in this list: " << std::endl;
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl;
  std::cout << "DO NOT PANIC - this is not a segfault" << std::endl;
}
