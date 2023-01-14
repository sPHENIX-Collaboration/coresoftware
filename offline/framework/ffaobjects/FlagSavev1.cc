#include "FlagSavev1.h"

#include <phool/PHFlag.h>

#include <cstdint>  // for uint64_t
#include <utility>  // for pair

class PHObject;

PHObject *
FlagSavev1::CloneMe() const
{
  FlagSavev1 *ret = new FlagSavev1();

  ret->intflag = this->intflag;
  ret->m_uint64flag_map = this->m_uint64flag_map;
  ret->doubleflag = this->doubleflag;
  ret->floatflag = this->floatflag;
  ret->stringflag = this->stringflag;

  return ret;
}

int FlagSavev1::isValid() const
{
  if (intflag.empty() &&
      doubleflag.empty() &&
      floatflag.empty() &&
      stringflag.empty() &&
      m_uint64flag_map.empty())
  {
    return 0;
  }
  return 1;
}

void FlagSavev1::identify(std::ostream &out) const
{
  out << "identify yourself: I am an FlagSavev1 Object" << std::endl;
  PrintIntFlag(out);
  Printuint64Flag(out);
  PrintDoubleFlag(out);
  PrintFloatFlag(out);
  PrintStringFlag(out);
  return;
}

int FlagSavev1::FillFromPHFlag(const PHFlag *flags, const bool clearold)
{
  if (clearold)
  {
    ClearAll();
  }
  int iret = FillIntFromPHFlag(flags);
  iret += Filluint64FromPHFlag(flags);
  iret += FillDoubleFromPHFlag(flags);
  iret += FillFloatFromPHFlag(flags);
  iret += FillStringFromPHFlag(flags);
  return iret;
}

int FlagSavev1::PutFlagsBack(PHFlag *flags, const bool overwrite)
{
  int iret = PutIntToPHFlag(flags, overwrite);
  iret += Putuint64ToPHFlag(flags, overwrite);
  iret += PutDoubleToPHFlag(flags, overwrite);
  iret += PutFloatToPHFlag(flags, overwrite);
  iret += PutStringToPHFlag(flags, overwrite);
  return iret;
}

int FlagSavev1::FillIntFromPHFlag(const PHFlag *flags)
{
  std::map<std::string, int>::const_iterator iter;
  const std::map<std::string, int> *intm = flags->IntMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
  {
    intflag[iter->first] = iter->second;
  }
  return 0;
}

int FlagSavev1::Filluint64FromPHFlag(const PHFlag *flags)
{
  std::map<std::string, uint64_t>::const_iterator iter;
  const std::map<std::string, uint64_t> *intm = flags->uint64Map();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
  {
    m_uint64flag_map[iter->first] = iter->second;
  }
  return 0;
}

int FlagSavev1::FillDoubleFromPHFlag(const PHFlag *flags)
{
  std::map<std::string, double>::const_iterator iter;
  const std::map<std::string, double> *intm = flags->DoubleMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
  {
    doubleflag[iter->first] = iter->second;
  }
  return 0;
}

int FlagSavev1::FillFloatFromPHFlag(const PHFlag *flags)
{
  std::map<std::string, float>::const_iterator iter;
  const std::map<std::string, float> *intm = flags->FloatMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
  {
    floatflag[iter->first] = iter->second;
  }
  return 0;
}

int FlagSavev1::FillStringFromPHFlag(const PHFlag *flags)
{
  std::map<std::string, std::string>::const_iterator iter;
  const std::map<std::string, std::string> *intm = flags->StringMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
  {
    std::string input(iter->second);
    stringflag[iter->first] = input;
  }
  return 0;
}
// general procedure, if overwrite is set, just set the flag
// if overwrite is false, call get flag with initialization which is used if it does not exist
int FlagSavev1::PutIntToPHFlag(PHFlag *flags, const bool overwrite)
{
  std::map<std::string, int>::const_iterator iter;
  for (iter = intflag.begin(); iter != intflag.end(); ++iter)
  {
    if (overwrite)
    {
      flags->set_IntFlag(iter->first, iter->second);
    }
    else
    {
      flags->get_IntFlag(iter->first, iter->second);
    }
  }
  return 0;
}

int FlagSavev1::Putuint64ToPHFlag(PHFlag *flags, const bool overwrite)
{
  std::map<std::string, uint64_t>::const_iterator iter;
  for (iter = m_uint64flag_map.begin(); iter != m_uint64flag_map.end(); ++iter)
  {
    if (overwrite)
    {
      flags->set_uint64Flag(iter->first, iter->second);
    }
    else
    {
      flags->get_uint64Flag(iter->first, iter->second);
    }
  }
  return 0;
}

int FlagSavev1::PutDoubleToPHFlag(PHFlag *flags, const bool overwrite)
{
  std::map<std::string, double>::const_iterator iter;
  for (iter = doubleflag.begin(); iter != doubleflag.end(); ++iter)
  {
    if (overwrite)
    {
      flags->set_DoubleFlag(iter->first, iter->second);
    }
    else
    {
      flags->get_DoubleFlag(iter->first, iter->second);
    }
  }
  return 0;
}

int FlagSavev1::PutFloatToPHFlag(PHFlag *flags, const bool overwrite)
{
  std::map<std::string, float>::const_iterator iter;
  for (iter = floatflag.begin(); iter != floatflag.end(); ++iter)
  {
    if (overwrite)
    {
      flags->set_FloatFlag(iter->first, iter->second);
    }
    else
    {
      flags->get_FloatFlag(iter->first, iter->second);
    }
  }
  return 0;
}

int FlagSavev1::PutStringToPHFlag(PHFlag *flags, const bool overwrite)
{
  std::map<std::string, std::string>::const_iterator iter;
  for (iter = stringflag.begin(); iter != stringflag.end(); ++iter)
  {
    if (overwrite)
    {
      flags->set_StringFlag(iter->first, iter->second);
    }
    else
    {
      flags->get_StringFlag(iter->first, iter->second);
    }
  }
  return 0;
}

void FlagSavev1::PrintIntFlag(std::ostream &os) const
{
  if (intflag.empty())
  {
    return;
  }
  std::map<std::string, int>::const_iterator iter;
  os << "Int Flags: " << std::endl;
  for (iter = intflag.begin(); iter != intflag.end(); ++iter)
  {
    os << iter->first << ": " << iter->second << std::endl;
  }
  return;
}

void FlagSavev1::Printuint64Flag(std::ostream &os) const
{
  if (m_uint64flag_map.empty())
  {
    return;
  }
  std::map<std::string, uint64_t>::const_iterator iter;
  os << "UInt64 Flags: " << std::endl;
  for (iter = m_uint64flag_map.begin(); iter != m_uint64flag_map.end(); ++iter)
  {
    os << iter->first << ": " << iter->second << std::endl;
  }
  return;
}

void FlagSavev1::PrintDoubleFlag(std::ostream &os) const
{
  if (doubleflag.empty())
  {
    return;
  }
  std::map<std::string, double>::const_iterator iter;
  os << "Double Flags: " << std::endl;
  for (iter = doubleflag.begin(); iter != doubleflag.end(); ++iter)
  {
    os << iter->first << ": " << iter->second << std::endl;
  }
  return;
}

void FlagSavev1::PrintFloatFlag(std::ostream &os) const
{
  if (floatflag.empty())
  {
    return;
  }
  std::map<std::string, float>::const_iterator iter;
  os << "Float Flags: " << std::endl;
  for (iter = floatflag.begin(); iter != floatflag.end(); ++iter)
  {
    os << iter->first << ": " << iter->second << std::endl;
  }
  return;
}

void FlagSavev1::PrintStringFlag(std::ostream &os) const
{
  if (stringflag.empty())
  {
    return;
  }
  std::map<std::string, std::string>::const_iterator iter;
  os << "String Flags: " << std::endl;
  for (iter = stringflag.begin(); iter != stringflag.end(); ++iter)
  {
    os << iter->first << ": " << iter->second << std::endl;
  }
  return;
}

void FlagSavev1::ClearAll()
{
  intflag.clear();
  doubleflag.clear();
  floatflag.clear();
  stringflag.clear();
  m_uint64flag_map.clear();
}
