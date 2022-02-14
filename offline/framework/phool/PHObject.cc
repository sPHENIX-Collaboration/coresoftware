#include "PHObject.h"

#include "phool.h"

#include <TSystem.h>

#include <iostream>

class TObject;

PHObject*
PHObject::CloneMe() const
{
  std::cout << PHWHERE << " CloneMe() not implemented by daugther class"
            << std::endl;
  return nullptr;
}

PHObject*
PHObject::clone() const
{
  std::cout << PHWHERE << " clone() is obsolete" << std::endl;
  return nullptr;
}

void PHObject::identify(std::ostream& os) const
{
  os << "identify yourself: I am a PHObject object" << std::endl;
}

void PHObject::Reset()
{
  // This is bad, previous events are not cleared
  std::cout << PHWHERE
            << " Reset() not implemented for " << this->GetName()
            << " PREVIOUS EVENT IS NOT CLEARED"
            << std::endl;
  std::cout << "You most probably miss a library, if so there are one or more messages at startup like:"
            << std::endl
            << std::endl;
  std::cout << "Warning in <TClass::TClass>: no dictionary for class XXX is available"
            << std::endl
            << std::endl;
  std::cout << "load the libraries which contain these classes and try again"
            << std::endl;
  gSystem->Exit(1);
}

int PHObject::isValid() const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid() not implemented by daughter class"
            << std::endl;
  return 0;
}



void PHObject::CopyFrom(const PHObject */*obj*/)
{
  std::cout << PHWHERE
            << " CopyFrom(const PHObject *obj) is not implemented" << std::endl;
  gSystem->Exit(1);
}

PHObject *PHObject::Clone(const char */*newname*/) const
{
std::cout << PHWHERE
	  << "You are overriding the TObject::Clone method which is not supported"  << std::endl;
  gSystem->Exit(1);
  return nullptr;
}

void PHObject::Copy (TObject &/*object*/) const
{
std::cout << PHWHERE
	  << "You are overriding the TObject::Copy method which is not supported"  << std::endl;
  gSystem->Exit(1);
}
