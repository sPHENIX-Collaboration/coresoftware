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
PHFlag::get_StringFlag(const std::string &name) const
{
  std::map<std::string, std::string>::const_iterator iter = m_StringFlagMap.find(name);
  if (iter != m_StringFlagMap.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getString: ERROR Unknown character Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return "";
}

const std::string
PHFlag::get_StringFlag(const std::string &name, const std::string &defaultval)
{
  std::map<std::string, std::string>::const_iterator iter = m_StringFlagMap.find(name);
  if (iter != m_StringFlagMap.end())
  {
    return iter->second;
  }
  else
  {
    set_StringFlag(name, defaultval);
    return get_StringFlag(name);
  }
}

void PHFlag::set_StringFlag(const std::string &name, const std::string &charstr)
{
  m_StringFlagMap[name] = charstr;
  return;
}

double PHFlag::get_DoubleFlag(const std::string &name) const
{
  std::map<std::string, double>::const_iterator iter = m_DoubleFlagMap.find(name);
  if (iter != m_DoubleFlagMap.end())
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
  std::map<std::string, double>::const_iterator iter = m_DoubleFlagMap.find(name);
  if (iter != m_DoubleFlagMap.end())
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
  m_DoubleFlagMap[name] = iflag;
  return;
}

float PHFlag::get_FloatFlag(const std::string &name) const
{
  std::map<std::string, float>::const_iterator iter = m_FloatFlagMap.find(name);
  if (iter != m_FloatFlagMap.end())
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
  std::map<std::string, float>::const_iterator iter = m_FloatFlagMap.find(name);
  if (iter != m_FloatFlagMap.end())
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
  m_FloatFlagMap[name] = iflag;
  return;
}

int PHFlag::get_IntFlag(const std::string &name) const
{
  std::map<std::string, int>::const_iterator iter = m_IntFlagMap.find(name);
  if (iter != m_IntFlagMap.end())
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
  std::map<std::string, int>::const_iterator iter = m_IntFlagMap.find(name);
  if (iter != m_IntFlagMap.end())
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
  m_IntFlagMap[name] = iflag;
  return;
}

uint64_t PHFlag::get_uint64Flag(const std::string &name) const
{
  std::map<std::string, uint64_t>::const_iterator iter = m_UInt64FlagMap.find(name);
  if (iter != m_UInt64FlagMap.end())
  {
    return iter->second;
  }
  std::cout << "PHFlag::getFlag: ERROR Unknown uint64 Flag " << name << std::endl;
  PrintStackTrace();
  std::cout << "The following flags are implemented: " << std::endl;
  Print();
  return 0;
}

uint64_t PHFlag::get_uint64Flag(const std::string &name, uint64_t defaultval)
{
  std::map<std::string, uint64_t>::const_iterator iter = m_UInt64FlagMap.find(name);
  if (iter != m_UInt64FlagMap.end())
  {
    return iter->second;
  }
  else
  {
    set_uint64Flag(name, defaultval);
    return get_uint64Flag(name);
  }
}

void PHFlag::set_uint64Flag(const std::string &name, const uint64_t iflag)
{
  m_UInt64FlagMap[name] = iflag;
  return;
}

void PHFlag::Print() const
{
  PrintIntFlags();
  Printuint64Flags();
  PrintFloatFlags();
  PrintDoubleFlags();
  PrintStringFlags();
  return;
}

void PHFlag::PrintIntFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
            << "Integer Flags:" << std::endl;
  std::map<std::string, int>::const_iterator intiter;
  for (intiter = m_IntFlagMap.begin(); intiter != m_IntFlagMap.end(); ++intiter)
  {
    std::cout << intiter->first << " is " << intiter->second << std::endl;
  }
  return;
}

void PHFlag::Printuint64Flags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
            << "uint64 Flags:" << std::endl;
  std::map<std::string, uint64_t>::const_iterator intiter;
  for (intiter = m_UInt64FlagMap.begin(); intiter != m_UInt64FlagMap.end(); ++intiter)
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
  for (doubleiter = m_DoubleFlagMap.begin(); doubleiter != m_DoubleFlagMap.end(); ++doubleiter)
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
  for (floatiter = m_FloatFlagMap.begin(); floatiter != m_FloatFlagMap.end(); ++floatiter)
  {
    std::cout << floatiter->first << " is " << floatiter->second << std::endl;
  }
  return;
}

void PHFlag::PrintStringFlags() const
{
  // loop over the map and print out the content (name and location in memory)
  std::cout << std::endl
            << "String Flags:" << std::endl;
  std::map<std::string, std::string>::const_iterator chariter;
  for (chariter = m_StringFlagMap.begin(); chariter != m_StringFlagMap.end(); ++chariter)
  {
    std::cout << chariter->first << " is " << chariter->second << std::endl;
  }
  return;
}

int PHFlag::FlagExist(const std::string &name) const
{
  std::map<std::string, int>::const_iterator iter = m_IntFlagMap.find(name);
  if (iter != m_IntFlagMap.end())
  {
    return 1;
  }
  std::map<std::string, uint64_t>::const_iterator uiter = m_UInt64FlagMap.find(name);
  if (uiter != m_UInt64FlagMap.end())
  {
    return 1;
  }
  std::map<std::string, float>::const_iterator fiter = m_FloatFlagMap.find(name);
  if (fiter != m_FloatFlagMap.end())
  {
    return 1;
  }
  std::map<std::string, double>::const_iterator diter = m_DoubleFlagMap.find(name);
  if (diter != m_DoubleFlagMap.end())
  {
    return 1;
  }
  std::map<std::string, std::string>::const_iterator citer = m_StringFlagMap.find(name);
  if (citer != m_StringFlagMap.end())
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
  uint64_t uivalue = 0;
  int uivaluecount = 0;
  std::string cvalue;
  int cvaluecount = 0;
  std::string junk;
  int junkcount = 0;

  std::ifstream infile(name);
  while (infile >> label)
  {
    std::cout << "Label " << label;
    if (label.substr(0, 1) == "S")
    {
      infile >> cvalue;
      cvaluecount++;
      set_StringFlag(label.substr(1, label.size() - 1), cvalue);
      std::cout << " type S read " << cvalue << std::endl;
    }
    else if (label.substr(0, 1) == "F")
    {
      infile >> fvalue;
      fvaluecount++;
      set_FloatFlag(label.substr(1, label.size() - 1), fvalue);
      std::cout << " type F read " << fvalue << std::endl;
    }
    else if (label.substr(0, 1) == "D")
    {
      infile >> dvalue;
      dvaluecount++;
      set_DoubleFlag(label.substr(1, label.size() - 1), dvalue);
      std::cout << " type D read " << dvalue << std::endl;
    }
    else if (label.substr(0, 1) == "I")
    {
      infile >> ivalue;
      ivaluecount++;
      set_IntFlag(label.substr(1, label.size() - 1), ivalue);
      std::cout << " type I read " << ivalue << std::endl;
    }
    else if (label.substr(0, 1) == "U")
    {
      infile >> uivalue;
      uivaluecount++;
      set_uint64Flag(label.substr(1, label.size() - 1), uivalue);
      std::cout << " type U read " << uivalue << std::endl;
    }
    else
    {
      infile >> junk;
      junkcount++;
      std::cout << " Junk read " << junk << std::endl;
    }
  }

  std::cout << "Read StringFlags(" << cvaluecount
            << ") FloatFlags(" << fvaluecount
            << ") DoubleFlags(" << dvaluecount
            << ") IntFlags(" << ivaluecount
            << ") uint64Flags(" << uivaluecount
            << ") JunkEntries(" << junkcount
            << ") from file " << name << std::endl;

  infile.close();
}

void PHFlag::WriteToFile(const std::string &name)
{
  std::ofstream outFile(name);
  // loop over the map and write out the content
  std::map<std::string, int>::const_iterator intiter;
  for (intiter = m_IntFlagMap.begin(); intiter != m_IntFlagMap.end(); ++intiter)
  {
    outFile << "I" << intiter->first << "\t" << intiter->second << std::endl;
  }

  std::map<std::string, uint64_t>::const_iterator uintiter;
  for (uintiter = m_UInt64FlagMap.begin(); uintiter != m_UInt64FlagMap.end(); ++uintiter)
  {
    outFile << "U" << uintiter->first << "\t" << uintiter->second << std::endl;
  }

  std::map<std::string, float>::const_iterator floatiter;
  for (floatiter = m_FloatFlagMap.begin(); floatiter != m_FloatFlagMap.end(); ++floatiter)
  {
    outFile << "F" << floatiter->first << "\t" << floatiter->second << std::endl;
  }

  int oldprecision = outFile.precision(15);
  std::map<std::string, double>::const_iterator doubleiter;
  for (doubleiter = m_DoubleFlagMap.begin(); doubleiter != m_DoubleFlagMap.end(); ++doubleiter)
  {
    outFile << "D" << doubleiter->first << "\t" << doubleiter->second << std::endl;
  }
  outFile.precision(oldprecision);

  std::map<std::string, std::string>::const_iterator chariter;
  for (chariter = m_StringFlagMap.begin(); chariter != m_StringFlagMap.end(); ++chariter)
  {
    outFile << "S" << chariter->first << "\t" << chariter->second << std::endl;
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

void PHFlag::ClearFlag(const std::string &name)
{
  std::map<std::string, int>::iterator iter = m_IntFlagMap.find(name);
  if (iter != m_IntFlagMap.end())
  {
    m_IntFlagMap.erase(iter);
  }
  std::map<std::string, uint64_t>::iterator uiter = m_UInt64FlagMap.find(name);
  if (uiter != m_UInt64FlagMap.end())
  {
    m_UInt64FlagMap.erase(uiter);
  }
  std::map<std::string, double>::iterator diter = m_DoubleFlagMap.find(name);
  if (diter != m_DoubleFlagMap.end())
  {
    m_DoubleFlagMap.erase(diter);
  }
  std::map<std::string, float>::iterator fiter = m_FloatFlagMap.find(name);
  if (fiter != m_FloatFlagMap.end())
  {
    m_FloatFlagMap.erase(fiter);
  }
  std::map<std::string, std::string>::iterator citer = m_StringFlagMap.find(name);
  if (citer != m_StringFlagMap.end())
  {
    m_StringFlagMap.erase(citer);
  }
}

void PHFlag::ClearAll()
{
  m_UInt64FlagMap.clear();
  m_IntFlagMap.clear();
  m_DoubleFlagMap.clear();
  m_FloatFlagMap.clear();
  m_StringFlagMap.clear();
  return;
}
