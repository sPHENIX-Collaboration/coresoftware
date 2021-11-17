#include "RunHeader.h"

#include <phool/phool.h>

#include <iostream>

class PHObject;

static int nowarning = 0;

PHObject*
RunHeader::CloneMe() const
{
  std::cout << "RunHeader::CloneMe() is not implemented in daugther class" << std::endl;
  return nullptr;
}

void RunHeader::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void RunHeader::identify(std::ostream& os) const
{
  os << "identify yourself: virtual RunHeader Object" << std::endl;
  return;
}

int RunHeader::isValid() const
{
  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
  return 0;
}

int RunHeader::get_RunNumber() const
{
  warning("get_RunNumber()");
  return -9999;
}

void RunHeader::set_RunNumber(const int /*run*/)
{
  warning("set_RunNumber(const int run)");
  return;
}

//
void RunHeader::NoWarning(const int i)
{
  if (i > 0)
  {
    std::cout << "RunHeader: switching virtual warnings OFF" << std::endl;
    nowarning = i;
  }
  else
  {
    std::cout << "RunHeader: switching virtual warnings ON" << std::endl;
    nowarning = 0;
  }
  return;
}

void RunHeader::warning(const std::string& funcname) const
{
  if (!nowarning)
  {
    std::cout << "Using virtual function RunHeader::" << funcname << std::endl;
  }
  return;
}
