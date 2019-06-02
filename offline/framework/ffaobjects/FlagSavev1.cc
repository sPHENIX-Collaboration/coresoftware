#include "FlagSavev1.h"

#include <phool/PHFlag.h>

#include <utility>         // for pair

using namespace std;

FlagSavev1 * 
FlagSavev1::clone() const
{
  FlagSavev1 * ret = new FlagSavev1();

  ret->intflag = this->intflag;
  ret->doubleflag = this->doubleflag;
  ret->floatflag = this->floatflag;
  ret->stringflag = this->stringflag;

  return ret;
}

int
FlagSavev1::isValid() const
{
  if (intflag.empty() &&
      doubleflag.empty() &&
      floatflag.empty() &&
      stringflag.empty())
    {
      return 0;
    }
  return 1;
}

void 
FlagSavev1::identify(ostream& out) const
{
  out << "identify yourself: I am an FlagSavev1 Object" << endl;
  PrintIntFlag(out) ;
  PrintDoubleFlag(out);
  PrintFloatFlag(out);
  PrintStringFlag(out);
  return ;
}

int
FlagSavev1::FillFromPHFlag(const PHFlag *flags)
{
  int iret = FillIntFromPHFlag(flags);
  iret += FillDoubleFromPHFlag(flags);
  iret += FillFloatFromPHFlag(flags);
  iret += FillCharFromPHFlag(flags);
  return iret;
}

int
FlagSavev1::PutFlagsBack(PHFlag *flags)
{
  int iret = PutIntToPHFlag(flags);
  iret += PutDoubleToPHFlag(flags);
  iret += PutFloatToPHFlag(flags);
  iret += PutCharToPHFlag(flags);
 return iret;
}


int
FlagSavev1::FillIntFromPHFlag(const PHFlag *flags)
{
  map<string, int>::const_iterator iter;
  const map<string, int> *intm = flags->IntMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
    {
      intflag[iter->first] = iter->second;
    }
  return 0;
}

int
FlagSavev1::FillDoubleFromPHFlag(const PHFlag *flags)
{
  map<string, double>::const_iterator iter;
  const map<string, double> *intm = flags->DoubleMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
    {
      doubleflag[iter->first] = iter->second;
    }
  return 0;
}

int
FlagSavev1::FillFloatFromPHFlag(const PHFlag *flags)
{
  map<string, float>::const_iterator iter;
  const map<string, float> *intm = flags->FloatMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
    {
      floatflag[iter->first] = iter->second;
    }
  return 0;
}

int
FlagSavev1::FillCharFromPHFlag(const PHFlag *flags)
{
  map<string, string>::const_iterator iter;
  const map<string, string> *intm = flags->CharMap();
  for (iter = intm->begin(); iter != intm->end(); ++iter)
    {
      string input(iter->second);
      stringflag[iter->first] = input;
    }
  return 0;
}

int
FlagSavev1::PutIntToPHFlag(PHFlag *flags)
{
  map<string, int>::const_iterator iter;
  for (iter = intflag.begin(); iter != intflag.end(); ++iter)
    {
      flags->set_IntFlag(iter->first,iter->second);
    }
  return 0;
}

int
FlagSavev1::PutDoubleToPHFlag(PHFlag *flags)
{
  map<string, double>::const_iterator iter;
  for (iter = doubleflag.begin(); iter != doubleflag.end(); ++iter)
    {
      flags->set_DoubleFlag(iter->first,iter->second);
    }
  return 0;
}

int
FlagSavev1::PutFloatToPHFlag(PHFlag *flags)
{
  map<string, float>::const_iterator iter;
  for (iter = floatflag.begin(); iter != floatflag.end(); ++iter)
    {
      flags->set_FloatFlag(iter->first,iter->second);
    }
  return 0;
}

int
FlagSavev1::PutCharToPHFlag(PHFlag *flags)
{
  map<string, string>::const_iterator iter;
  for (iter = stringflag.begin(); iter != stringflag.end(); ++iter)
    {
      flags->set_CharFlag(iter->first,iter->second);
    }
  return 0;
}

void
FlagSavev1::PrintIntFlag(std::ostream& os) const
{
  if (intflag.empty())
    {
      return ;
    }
  map<string, int>::const_iterator iter;
  os << "Int Flags: " << endl;
  for (iter = intflag.begin(); iter != intflag.end(); ++iter)
    {
      os << iter->first << ": " << iter->second << endl;
    }
  return ;
}

void
FlagSavev1::PrintDoubleFlag(std::ostream& os) const
{
  if (doubleflag.empty())
    {
      return ;
    }
  map<string, double>::const_iterator iter;
  os << "Double Flags: " << endl;
  for (iter = doubleflag.begin(); iter != doubleflag.end(); ++iter)
    {
      os << iter->first << ": " << iter->second << endl;
    }
  return ;
}

void
FlagSavev1::PrintFloatFlag(std::ostream& os) const
{
  if (floatflag.empty())
    {
      return ;
    }
  map<string, float>::const_iterator iter;
  os << "Float Flags: " << endl;
  for (iter = floatflag.begin(); iter != floatflag.end(); ++iter)
    {
      os << iter->first << ": " << iter->second << endl;
    }
  return ;
}

void
FlagSavev1::PrintStringFlag(std::ostream& os) const
{
  if (stringflag.empty())
    {
      return ;
    }
  map<string, string>::const_iterator iter;
  os << "String Flags: " << endl;
    for (iter = stringflag.begin(); iter != stringflag.end(); ++iter)
      {
        os << iter->first << ": " << iter->second << endl;
      }
    return ;
  }
