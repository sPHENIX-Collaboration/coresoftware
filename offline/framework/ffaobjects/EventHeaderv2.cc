#include "EventHeaderv2.h"

#include <iostream>

void EventHeaderv2::Reset()
{ 
  EventHeaderv1::Reset();
  m_bunchCrossing = 0; 
}

void EventHeaderv2::identify(std::ostream& out) const
{
  out << "identify yourself: I am an EventHeaderv2 Object" << std::endl;
  const auto timestamp( get_TimeStamp() );
  out 
    << "Run Number: " << get_RunNumber() 
    << ", Event no: " << get_EvtSequence() 
    << ", Type: " << get_EvtType() 
    << ", DAQ arrival time: " << ctime(&timestamp)
    << ", bunch crossing: " << m_bunchCrossing
    << std::endl;
}

