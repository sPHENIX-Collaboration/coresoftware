#include "SyncObject.h"

#include <phool/phool.h>

#include <iostream>

using namespace std;

void SyncObject::Reset()
{
  cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << endl;
  return ;
}

void SyncObject::identify(ostream& os) const
{
  os << "identify yourself: virtual SyncObject Object" << endl;
  return ;
}

int SyncObject::isValid() const
{
  cout << PHWHERE << "isValid not implemented by daughter class" << endl;
  return 0;
}

SyncObject*
SyncObject::clone() const
{
  cout << "SyncObject::clone() not implemented by daughter class" << endl;
  return nullptr;
}

SyncObject&
SyncObject::operator=(const SyncObject& source)
{
  if ( this != &source )
    {
      EventCounter(source.EventCounter());
      EventNumber(source.EventNumber());
      RunNumber(source.RunNumber());
      SegmentNumber(source.SegmentNumber());
    }
  return *this;
}

int SyncObject::Different(const SyncObject *other) const
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
