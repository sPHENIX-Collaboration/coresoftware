#include "MbdRawHit.h"

#include <iostream>

void MbdRawHit::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdRawHit object" << std::endl;
}
