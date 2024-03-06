
#include "LL1ReturnCodes.h"
#include "LL1Out.h"

#include <cmath>
#include <iostream>

ClassImp(LL1Out)

LL1Out::LL1Out()
{
  Init();
}
LL1Out::~LL1Out()
{

}

//______________________________________
void LL1Out::Init()
{
  _trigger_type = "NONE";
}


//______________________________________
void LL1Out::Reset()
{
  Init();
}

//______________________________________
void LL1Out::identify(std::ostream& out)
{
  out << "identify yourself: I am a LL1Out object" << std::endl;
  out << "triggertype: " << _trigger_type << std::endl;
}

int LL1Out::isValid()
{

  return 1;
}
