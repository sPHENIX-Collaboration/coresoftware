#include "SyncObjectv2.h"
#include <iostream>


ClassImp(SyncObjectv2)

using namespace std;

SyncObjectv2::SyncObjectv2():
  segmentnumber(-999999)
{
  return;
}

SyncObjectv2::SyncObjectv2(const SyncObject& source)
{
  EventCounter(source.EventCounter());
  EventNumber(source.EventNumber());
  RunNumber(source.RunNumber());
  SegmentNumber(source.SegmentNumber());
}

void SyncObjectv2::Reset()
{
  segmentnumber = -999999;
  SyncObjectv1::Reset();
  return;
}

void SyncObjectv2::identify(ostream& out) const
{
  out << "identify yourself: I am an SyncObjectv2 Object" << endl;
  out << "Event no: " <<  eventnumber
      << ", Counter: " << eventcounter
      << ", Segment Number: " << segmentnumber
      << ", Run Number: " << runnumber
      << endl;

  return;
}

