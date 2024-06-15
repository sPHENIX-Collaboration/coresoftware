#include "TriggerPrimitiveContainerv1.h"
#include "TriggerDefs.h"

#include <phool/PHObject.h>

#include <TClonesArray.h>

#include <cassert>
#include <cmath>
#include <iostream>

TriggerPrimitiveContainerv1::TriggerPrimitiveContainerv1(const TriggerDefs::TriggerId tid, const TriggerDefs::DetectorId did)
{
  m_triggerid = tid;
  m_detectorid = did;

  // Make the primitives.
  int nprimitives = 0;
  
  int nsums = 0;
  

  if (tid == TriggerDefs::noneTId)
    {
      if (did == TriggerDefs::DetectorId::emcalDId)
	{
	  nprimitives = 384;
	}
      else if (did == TriggerDefs::DetectorId::hcalinDId || did == TriggerDefs::DetectorId::hcaloutDId || did == TriggerDefs::DetectorId::hcalDId)
	{
	  nprimitives = 24;
	}
      m_primitiveid = TriggerDefs::PrimitiveId::calPId;
      nsums = 16;
    }
  else if (tid == TriggerDefs::jetTId || tid == TriggerDefs::photonTId)
    {
      nprimitives = 16;

      m_primitiveid = TriggerDefs::PrimitiveId::jetPId;
      nsums = 24;
 
    }
  else if (tid == TriggerDefs::pairTId)
    {
      nprimitives = 16;
      m_primitiveid = TriggerDefs::PrimitiveId::pairPId;
      nsums = 16;
    }
  else if (tid == TriggerDefs::mbdTId)
    {
      nprimitives = 4;
      m_primitiveid = TriggerDefs::PrimitiveId::mbdPId;
      nsums = 13;
    }
  // add a primitive for all primitives
  for (int ip = 0 ;ip < nprimitives; ip++)
    {
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(m_triggerid, m_detectorid, m_primitiveid, ip);
      TriggerPrimitive *primitive = new TriggerPrimitivev1(primkey);
      for (int is = 0; is < nsums; is++)
	{ 
	  TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(m_triggerid, m_detectorid, m_primitiveid, ip, is);
	  std::vector<unsigned int> *t_sum = new std::vector<unsigned int>();
	  primitive->add_sum(sumkey, t_sum);
	}
      add_primitive(primkey, primitive);
    }
  return;
}

TriggerPrimitiveContainerv1::~TriggerPrimitiveContainerv1()
= default;

//______________________________________
void TriggerPrimitiveContainerv1::Reset()
{
  for ( const auto &primitive : _primitives)
  {
    primitive.second->Reset();
  }
}

//______________________________________
void TriggerPrimitiveContainerv1::identify(std::ostream& out) const
{
  out << " I am a Trigger Primitive Container :: Trigger id: " << std::hex << m_triggerid << std::endl;

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
