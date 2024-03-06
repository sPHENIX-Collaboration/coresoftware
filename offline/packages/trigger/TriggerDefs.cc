#include "TriggerDefs.h"
#include <bitset>
#include <cstring>
#include <calobase/TowerInfoDefs.h>

uint32_t
TriggerDefs::getTriggerKey(const TriggerDefs::TriggerId triggerId)
{
  uint32_t tmp = triggerId << kBitShiftTriggerId;
  return tmp;
}

uint32_t
TriggerDefs::getTriggerPrimKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid)
{
  uint32_t tmp = triggerId;
  uint32_t key = triggerId << kBitShiftTriggerId;
  tmp = detectorId;
  key |= tmp << kBitShiftDetectorId;
  tmp = primitiveId;
  key |= tmp << kBitShiftPrimitiveId;
  tmp = primlocid;
  key |= tmp << kBitShiftPrimitiveLocId;
  return key;
}

uint32_t
TriggerDefs::getTriggerSumKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId, const TriggerDefs::PrimitiveId primitiveId, const uint16_t primlocid, const uint16_t sumlocid)
{
  uint32_t tmp = triggerId;
  uint32_t key = triggerId << kBitShiftTriggerId;
  tmp = detectorId;
  key |= tmp << kBitShiftDetectorId;
  tmp = primitiveId;
  key |= tmp << kBitShiftPrimitiveId;
  tmp = primlocid;
  key |= tmp << kBitShiftPrimitiveLocId;
  tmp = sumlocid;
  key |= tmp << kBitShiftSumLocId;
  return key;
}

uint32_t
TriggerDefs::getTriggerId_from_TriggerKey(const TriggerDefs::TriggerKey triggerkey)
{
  uint32_t tmp = (triggerkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t
TriggerDefs::getTriggerId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t 
TriggerDefs::getTriggerId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftTriggerId);
  return tmp;
}

uint32_t
TriggerDefs::getDetectorId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftDetectorId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getDetectorId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftDetectorId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getPrimitiveId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp =  (triggerprimkey >> kBitShiftPrimitiveId) & 0xf;
  return tmp;
}

uint32_t 
TriggerDefs::getPrimitiveId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp =  (triggersumkey >> kBitShiftPrimitiveId) & 0xf;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint16_t tmp =  (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitiveLocId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp =  (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  return tmp;
}

uint16_t 
TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId= getDetectorId_from_TriggerPrimKey(triggerprimkey);
  uint32_t primId= getPrimitiveId_from_TriggerPrimKey(triggerprimkey);
  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp/3;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      if (primId == TriggerDefs::PrimitiveId::calPId)
	{
	  return tmp/12;
	  break;
	}
      else if (primId == TriggerDefs::PrimitiveId::jetPId)
	{
	  return tmp;
	  break;
	}
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);
  uint32_t primId= getPrimitiveId_from_TriggerPrimKey(triggersumkey);
  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp/3;
      break;

    case TriggerDefs::DetectorId::emcalDId :
      if (primId == TriggerDefs::PrimitiveId::calPId)
	{
	  return tmp/12;
	  break;
	}
      if (primId == TriggerDefs::PrimitiveId::jetPId)
	{
	  return tmp;
	  break;
	}

      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId= getDetectorId_from_TriggerPrimKey(triggerprimkey);
  uint32_t primId= getPrimitiveId_from_TriggerPrimKey(triggerprimkey);

  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp%3;
      break;
    case TriggerDefs::DetectorId::emcalDId :

      if (primId == TriggerDefs::PrimitiveId::calPId)
	{
	  return tmp%12;
	  break;
	}
      else if (primId == TriggerDefs::PrimitiveId::jetPId)
	{
	  return 0;
	  break;
	}
      
    default :
      return UINT16_MAX;
      break;
    }

  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);
  uint32_t primId= getPrimitiveId_from_TriggerSumKey(triggersumkey);
  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ff;
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp%3;
      break;

    case TriggerDefs::DetectorId::emcalDId :
      if (primId == TriggerDefs::PrimitiveId::calPId)
	{
	  return tmp%12;
	  break;
	}
      if (primId == TriggerDefs::PrimitiveId::jetPId)
	{
	  return 0;
	  break;
	}
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint16_t 
TriggerDefs::getSumLocId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp = (triggersumkey >> kBitShiftSumLocId) & 0x1f;
  return tmp;
}
uint16_t 
TriggerDefs::getSumPhiId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;


      return tmp/4;
      break;
    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp/4;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp/4;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}
uint16_t 
TriggerDefs::getSumEtaId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId= getDetectorId_from_TriggerSumKey(triggersumkey);

  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);
  
  switch (detId)
    {
    case TriggerDefs::DetectorId::mbdDId :
      return UINT16_MAX;
      break;

    case TriggerDefs::DetectorId::hcalinDId :
    case TriggerDefs::DetectorId::hcaloutDId :
    case TriggerDefs::DetectorId::hcalDId :
      return tmp%4;
      break;
    case TriggerDefs::DetectorId::emcalDId :
      return tmp%4;
      break;
    default :
      return UINT16_MAX;
      break;
    }
  return UINT16_MAX;
}

uint32_t 
TriggerDefs::GetTowerInfoKey( const TriggerDefs::DetectorId detId, const uint16_t iprim, const uint16_t isum, const uint16_t itower )
{
  unsigned int phibin = 0;
  unsigned int etabin = 0;
  switch (detId)
    {
    case TriggerDefs::DetectorId::emcalDId:
      phibin = 8*(iprim/12) + 2*(isum/4) + (itower/2);
      etabin = 8*(iprim%12) + 2*(isum%4) + (itower%2);
      if (phibin < 0 || phibin > 255 || etabin < 0 || etabin > 95) 
	{
	  std::cout << "PHIBIN eTABIN BADD: "<<phibin << " " <<etabin<<std::endl;
	  return UINT16_MAX; 
	  break;
	}
      return TowerInfoDefs::encode_emcal(etabin, phibin);	  ;
      break;

    default :
      return UINT16_MAX; 
      break;
    }

  return UINT16_MAX;
}
TriggerDefs::TriggerId TriggerDefs::GetTriggerId(std::string trigger)
  {

    if (strcmp(trigger.c_str(), "NONE") == 0) return TriggerDefs::TriggerId::noneTId;
    else if (strcmp(trigger.c_str(), "MBD") == 0) return TriggerDefs::TriggerId::mbdTId;
    else if (strcmp(trigger.c_str(), "JET") == 0) return TriggerDefs::TriggerId::jetTId;
    else if (strcmp(trigger.c_str(), "PHOTON") == 0) return TriggerDefs::TriggerId::photonTId;
    else if (strcmp(trigger.c_str(), "PAIR") == 0) return TriggerDefs::TriggerId::pairTId;
    else if (strcmp(trigger.c_str(), "COSMIC") == 0) return TriggerDefs::TriggerId::cosmicTId;
    else if (strcmp(trigger.c_str(), "COSMIC_COIN") == 0) return TriggerDefs::TriggerId::cosmic_coinTId;
    
    return TriggerDefs::TriggerId::noneTId;
    
  }

TriggerDefs::DetectorId TriggerDefs::GetDetectorId(std::string detector)
  {

    if (strcmp(detector.c_str(), "NONE") == 0) return TriggerDefs::DetectorId::noneDId;
    else if (strcmp(detector.c_str(), "MBD") == 0) return TriggerDefs::DetectorId::mbdDId;
    else if (strcmp(detector.c_str(), "HCALIN") == 0) return TriggerDefs::DetectorId::hcalinDId;
    else if (strcmp(detector.c_str(), "HCALOUT") == 0) return TriggerDefs::DetectorId::hcaloutDId;
    else if (strcmp(detector.c_str(), "EMCAL") == 0) return TriggerDefs::DetectorId::emcalDId;
    
    return TriggerDefs::DetectorId::noneDId;
    
  }
TriggerDefs::PrimitiveId TriggerDefs::GetPrimitiveId(std::string primitive)
{

    if (strcmp(primitive.c_str(), "NONE") == 0) return TriggerDefs::PrimitiveId::nonePId;
    else if (strcmp(primitive.c_str(), "MBD") == 0) return TriggerDefs::PrimitiveId::mbdPId;
    else if (strcmp(primitive.c_str(), "HCALIN") == 0) return TriggerDefs::PrimitiveId::calPId;
    else if (strcmp(primitive.c_str(), "HCALOUT") == 0) return TriggerDefs::PrimitiveId::calPId;
    else if (strcmp(primitive.c_str(), "HCAL") == 0) return TriggerDefs::PrimitiveId::calPId;
    else if (strcmp(primitive.c_str(), "EMCAL") == 0) return TriggerDefs::PrimitiveId::calPId;
    else if (strcmp(primitive.c_str(), "JET") == 0) return TriggerDefs::PrimitiveId::jetPId;
    
    return TriggerDefs::PrimitiveId::nonePId;

  }
