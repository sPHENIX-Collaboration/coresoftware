#include "RunHeader.h"

#include <phool/phool.h>

#include <iostream>

class PHObject;

using namespace std;

static int nowarning = 0;

PHObject *
RunHeader::CloneMe() const
{
  cout << "RunHeader::CloneMe() is not implemented in daugther class" << endl;
  return nullptr;
}

void
RunHeader::Reset()
{
  cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << endl;
  return ;
}

void
RunHeader::identify(ostream& os) const
{
  os << "identify yourself: virtual RunHeader Object" << endl;
  return ;
}

int
RunHeader::isValid() const
{
  cout << PHWHERE << "isValid not implemented by daughter class" << endl;
  return 0;
}

int
RunHeader::get_RunNumber() const
{
  warning("get_RunNumber()");
  return -9999;
}

void
RunHeader::set_RunNumber(const int /*run*/)
{
  warning("set_RunNumber(const int run)");
  return ;
}

double
RunHeader::get_Bfield() const
{
  warning("get_Bfield()");
  return -9999.9;
}

void
RunHeader::set_Bfield(const double /*rval*/)
{
  warning("set_Bfield(const double rval)");
  return;
}

time_t
RunHeader::get_TimeStart() const
{
  warning("get_TimeStart()");
  return 0;
}

void
RunHeader::set_TimeStart(const time_t /*ival*/)
{
  warning("set_TimeStart(const time_t ival)");
  return;
}

time_t
RunHeader::get_TimeStop() const
{
  warning("get_TimeStop()");
  return 0;
}

void
RunHeader::set_TimeStop(const time_t /*ival*/)
{
  warning("set_TimeStop(const time_t ival)");
  return;
}

// 
void
RunHeader::NoWarning(const int i)
{
  if (i > 0)
    {
      cout << "RunHeader: switching virtual warnings OFF" << endl;
      nowarning = i;
    }
  else
    {
      cout << "RunHeader: switching virtual warnings ON" << endl;
      nowarning = 0;
    }
  return ;
}

void
RunHeader::warning(const char *funcname) const
{
  if (! nowarning)
    {
      cout << "Using virtual function RunHeader::" << funcname << endl;
    }
  return ;
}
