
#include "TriggerPrimitiveContainer.h"

#include <cmath>
#include <iostream>

ClassImp(TriggerPrimitiveContainer)

TriggerPrimitiveContainer::TriggerPrimitiveContainer()
{
}
TriggerPrimitiveContainer::~TriggerPrimitiveContainer()
{

}

//______________________________________
void TriggerPrimitiveContainer::Reset()
{

}

void TriggerPrimitiveContainer::add_primitive(TriggerDefs::TriggerKey /*key*/, TriggerPrimitive* /*primitive*/)
{

}



//______________________________________
void TriggerPrimitiveContainer::identify(std::ostream& out)
{
  out << "identify yourself: I am a TriggerPrimitiveContainer object" << std::endl;

}

int TriggerPrimitiveContainer::isValid()
{
  return (!_primitives.empty());
}
