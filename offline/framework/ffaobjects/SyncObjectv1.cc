#include "SyncObjectv1.h"

SyncObjectv1::SyncObjectv1(const SyncObject& source)
{
  EventCounter(source.EventCounter());
  EventNumber(source.EventNumber());
  RunNumber(source.RunNumber());
  SegmentNumber(source.SegmentNumber());
}

void SyncObjectv1::Reset()
{
  eventnumber = 0;
  runnumber = 0;
  eventcounter = 0;
  segmentnumber = -999999;
  return;
}

void SyncObjectv1::identify(std::ostream& out) const
{
  out << "identify yourself: I am an SyncObjectv1 Object" << std::endl;
  out << "Event no: " << eventnumber
      << ", Counter: " << eventcounter
      << ", Segment Number: " << segmentnumber
      << ", Run Number: " << runnumber
      << std::endl;

  return;
}

int SyncObjectv1::isValid() const
{
  return ((eventcounter) ? 1 : 0);  // return 1 if eventcounter is not zero
}
