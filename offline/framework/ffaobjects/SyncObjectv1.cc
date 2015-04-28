#include "SyncObjectv1.h"

using namespace std;
ClassImp(SyncObjectv1)

SyncObjectv1::SyncObjectv1():
  eventcounter(0),
  eventnumber(-999999),
  runnumber(-999999)
{
  return;
}

SyncObjectv1::SyncObjectv1(const SyncObject& source)
{
  EventCounter(source.EventCounter());
  EventNumber(source.EventNumber());
  RunNumber(source.RunNumber());
}

void SyncObjectv1::Reset()
{
  eventnumber = -999999;
  runnumber = -999999;
  eventcounter = 0;
  return;
}

void SyncObjectv1::identify(ostream& out) const
{
  out << "identify yourself: I am an SyncObjectv1 Object" << endl;
  out << "Event no: " <<  eventnumber
      << ", Counter: " << eventcounter
      << ", Run Number: " << runnumber
      << endl;

  return;
}

int SyncObjectv1::isValid() const
{
  return((eventcounter) ? 1:0); // return 1 if eventcounter is not zero
}

