#include "EventHeader_Prototype2.h"

using namespace std;

ClassImp(EventHeader_Prototype2)

EventHeader_Prototype2::EventHeader_Prototype2()
{
  Reset();
  return;
}

void 
EventHeader_Prototype2::Reset()
{
  TimeStamp = 0;
  EvtSequence = -999999;
  EvtType = -999999;
  return;
}

void 
EventHeader_Prototype2::identify(ostream& out) const
{
  out << "identify yourself: I am an EventHeader_Prototype2 Object" << endl;
  out << "Event no: " << EvtSequence 
      << ", Type: " << EvtType 
      << ", RCDAQ arrival time: " << ctime(&TimeStamp)
      << endl;

  return;
}

int 
EventHeader_Prototype2::isValid() const
{
  return((TimeStamp) ? 1:0); // return 1 if TimeStamp is not zero
}
