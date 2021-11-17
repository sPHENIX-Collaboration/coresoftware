#include "SyncObject.h"

#include <phool/phool.h>

#include <iostream>

class PHObject;

void SyncObject::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void SyncObject::identify(std::ostream& os) const
{
  os << "identify yourself: virtual SyncObject Object" << std::endl;
  return;
}

int SyncObject::isValid() const
{
  std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
  return 0;
}

PHObject*
SyncObject::CloneMe() const
{
  std::cout << "SyncObject::CloneMe() not implemented by daughter class" << std::endl;
  return nullptr;
}

SyncObject&
SyncObject::operator=(const SyncObject& source)
{
  if (this != &source)
  {
    EventCounter(source.EventCounter());
    EventNumber(source.EventNumber());
    RunNumber(source.RunNumber());
    SegmentNumber(source.SegmentNumber());
  }
  return *this;
}

int SyncObject::Different(const SyncObject* other) const
{
  int iret = 0;
  if (EventNumber() != other->EventNumber())
  {
    iret += 0x1;
  }
  if (RunNumber() != other->RunNumber())
  {
    iret |= 0x2;
  }
  return iret;
}
