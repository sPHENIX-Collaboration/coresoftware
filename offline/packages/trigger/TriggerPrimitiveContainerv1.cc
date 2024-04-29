#include "TriggerPrimitiveContainerv1.h"
#include "TriggerDefs.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

#include <cassert>
#include <cmath>
#include <iostream>

TriggerPrimitiveContainerv1::TriggerPrimitiveContainerv1(const std::string& triggertype)
{
  m_triggerkey = TriggerDefs::getTriggerKey(TriggerDefs::GetTriggerId(triggertype));
}

TriggerPrimitiveContainerv1::~TriggerPrimitiveContainerv1()
{
  TriggerPrimitiveContainerv1::Reset();
}

//______________________________________
void TriggerPrimitiveContainerv1::Reset()
{
  while (_primitives.begin() != _primitives.end())
  {
    delete _primitives.begin()->second;
    _primitives.erase(_primitives.begin());
  }
}

//______________________________________
void TriggerPrimitiveContainerv1::identify(std::ostream& out) const
{
  out << " I am a Trigger Primitive Container :: Trigger key: " << std::hex << m_triggerkey << std::endl;

  out << " TriggerId: " << TriggerDefs::getTriggerId_from_TriggerKey(m_triggerkey) << std::endl;

  ConstRange range = getTriggerPrimitives();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    (*iter).second->identify();
  }
}

int TriggerPrimitiveContainerv1::isValid() const
{
  return (!_primitives.empty());
}

TriggerPrimitive* TriggerPrimitiveContainerv1::get_primitive_at_key(TriggerDefs::TriggerPrimKey key)
{
  TriggerDefs::TriggerPrimKey pkey = (key & 0xffffffe0U);
  if (!_primitives[pkey])
  {
    return nullptr;
  }

  return _primitives[pkey];
}

void TriggerPrimitiveContainerv1::add_primitive(TriggerDefs::TriggerPrimKey key, TriggerPrimitive* prim)
{
  _primitives[key] = prim;
}
TriggerPrimitiveContainerv1::ConstRange TriggerPrimitiveContainerv1::getTriggerPrimitives() const
{
  return make_pair(_primitives.begin(), _primitives.end());
}

TriggerPrimitiveContainerv1::Range TriggerPrimitiveContainerv1::getTriggerPrimitives()
{
  return make_pair(_primitives.begin(), _primitives.end());
}
