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
  out
      << "Run Number: " << get_RunNumber()
      << ", Event no: " << get_EvtSequence()
      << ", Type: " << get_EvtType()
      << ", bunch crossing: " << m_bunchCrossing
      << std::endl;
}
