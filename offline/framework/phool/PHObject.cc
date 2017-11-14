#include "PHObject.h"
#include <cstdlib>
#include <iostream>
#include "phool.h"

PHObject::PHObject()
  : split(99)
  , bufSize(32000)
{
  return;
}

PHObject*
PHObject::clone() const
{
  std::cout << PHWHERE << " clone() not implemented by daugther class"
            << std::endl;
  return 0;
}

void PHObject::identify(std::ostream& out) const
{
  out << "identify yourself: I am a PHObject object" << std::endl;
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
  exit(1);
}

int PHObject::isValid() const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid() not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isValid(const float) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid(const float f) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isValid(const double) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid(const double) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isValid(const int) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid(const int) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isValid(const unsigned int) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isValid(const unsigned int) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isImplemented(const float) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isImplemented(const float f) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isImplemented(const double) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isImplemented(const double) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isImplemented(const int) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isImplemented(const int) not implemented by daughter class"
            << std::endl;
  return 0;
}

int PHObject::isImplemented(const unsigned int) const
{
  // give warning if this method is not implemented
  std::cout << PHWHERE
            << " isImplemented(const unsigned int) not implemented by daughter class"
            << std::endl;
  return 0;
}

void PHObject::CopyContent(PHObject* obj)
{
  std::cout << PHWHERE
            << " CopyContent(PHObject *obj) is not implemented" << std::endl;
  exit(1);
}
