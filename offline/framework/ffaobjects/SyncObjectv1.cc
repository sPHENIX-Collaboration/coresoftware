#include "SyncObjectv1.h"

using namespace std;

SyncObjectv1::SyncObjectv1():
  eventcounter(0),
  eventnumber(-999999),
  runnumber(-999999),
  segmentnumber(-999999)
{
  return;
}

SyncObjectv1::SyncObjectv1(const SyncObject& source)
{
  EventCounter(source.EventCounter());
  EventNumber(source.EventNumber());
  RunNumber(source.RunNumber());
  SegmentNumber(source.SegmentNumber());
}

void SyncObjectv1::Reset()
{
  eventnumber = -999999;
  runnumber = -999999;
  eventcounter = 0;
  segmentnumber = -999999;
  return;
}

void SyncObjectv1::identify(ostream& out) const
{
  out << "identify yourself: I am an SyncObjectv1 Object" << endl;
  out << "Event no: " <<  eventnumber
      << ", Counter: " << eventcounter
      << ", Segment Number: " << segmentnumber
      << ", Run Number: " << runnumber
      << endl;

  return;
}

int SyncObjectv1::isValid() const
{
  return((eventcounter) ? 1:0); // return 1 if eventcounter is not zero
}

