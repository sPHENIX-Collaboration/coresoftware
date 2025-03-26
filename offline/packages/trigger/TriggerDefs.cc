#include "TriggerDefs.h"

#include <calobase/TowerInfoDefs.h>

uint32_t
TriggerDefs::getTriggerKey(const TriggerDefs::TriggerId triggerId)
{
  uint32_t tmp = triggerId << kBitShiftTriggerId;
  return tmp;
}

uint32_t
TriggerDefs::getTriggerKey(const TriggerDefs::TriggerId triggerId, const TriggerDefs::DetectorId detectorId)
{
  uint32_t tmp = triggerId << kBitShiftTriggerId;
  tmp |= (unsigned int) detectorId << kBitShiftDetectorId;
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

TriggerDefs::TriggerId
TriggerDefs::getTriggerId_from_TriggerKey(const TriggerDefs::TriggerKey triggerkey)
{
  uint32_t tmp = (triggerkey >> kBitShiftTriggerId);
  return static_cast<TriggerDefs::TriggerId>(tmp);
}

TriggerDefs::TriggerId
TriggerDefs::getTriggerId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp = (triggerprimkey >> kBitShiftTriggerId);
  return static_cast<TriggerDefs::TriggerId>(tmp);
}

TriggerDefs::TriggerId
TriggerDefs::getTriggerId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp = (triggersumkey >> kBitShiftTriggerId);
  return static_cast<TriggerDefs::TriggerId>(tmp);
}

TriggerDefs::DetectorId
TriggerDefs::getDetectorId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp = (triggerprimkey >> kBitShiftDetectorId) & 0xfU;
  return static_cast<TriggerDefs::DetectorId>(tmp);
}

TriggerDefs::DetectorId
TriggerDefs::getDetectorId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp = (triggersumkey >> kBitShiftDetectorId) & 0xfU;
  return static_cast<TriggerDefs::DetectorId>(tmp);
}

TriggerDefs::PrimitiveId
TriggerDefs::getPrimitiveId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t tmp = (triggerprimkey >> kBitShiftPrimitiveId) & 0xfU;
  return static_cast<TriggerDefs::PrimitiveId>(tmp);
}

TriggerDefs::PrimitiveId
TriggerDefs::getPrimitiveId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t tmp = (triggersumkey >> kBitShiftPrimitiveId) & 0xfU;
  return static_cast<TriggerDefs::PrimitiveId>(tmp);
}

uint16_t
TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ffU;
  return tmp;
}

uint16_t
TriggerDefs::getPrimitiveLocId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ffU;
  return tmp;
}

uint16_t
TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId = getDetectorId_from_TriggerPrimKey(triggerprimkey);
  uint32_t primId = getPrimitiveId_from_TriggerPrimKey(triggerprimkey);
  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ffU;

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    if (detId == TriggerDefs::DetectorId::hcalinDId ||
        detId == TriggerDefs::DetectorId::hcaloutDId ||
        detId == TriggerDefs::DetectorId::hcalDId)
    {
      return tmp / 3;
    }
    if (detId == TriggerDefs::DetectorId::emcalDId)
    {
      return tmp / 12;
    }
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return tmp;
  }

  return UINT16_MAX;
}

uint16_t
TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId = getDetectorId_from_TriggerSumKey(triggersumkey);
  uint32_t primId = getPrimitiveId_from_TriggerSumKey(triggersumkey);
  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ffU;

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    if (detId == TriggerDefs::DetectorId::hcalinDId ||
        detId == TriggerDefs::DetectorId::hcaloutDId ||
        detId == TriggerDefs::DetectorId::hcalDId)
    {
      return tmp / 3;
    }
    if (detId == TriggerDefs::DetectorId::emcalDId)
    {
      return tmp / 12;
    }
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return tmp;
  }
  return UINT16_MAX;
}

uint16_t
TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(const TriggerDefs::TriggerPrimKey triggerprimkey)
{
  uint32_t detId = getDetectorId_from_TriggerPrimKey(triggerprimkey);
  uint32_t primId = getPrimitiveId_from_TriggerPrimKey(triggerprimkey);

  uint16_t tmp = (triggerprimkey >> kBitShiftPrimitiveLocId) & 0x1ffU;

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    if (detId == TriggerDefs::DetectorId::hcalinDId ||
        detId == TriggerDefs::DetectorId::hcaloutDId ||
        detId == TriggerDefs::DetectorId::hcalDId)
    {
      return tmp % 3;
    }
    if (detId == TriggerDefs::DetectorId::emcalDId)
    {
      return tmp % 12;
    }
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return 0;
  }

  return UINT16_MAX;
}

uint16_t
TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t detId = getDetectorId_from_TriggerSumKey(triggersumkey);
  uint32_t primId = getPrimitiveId_from_TriggerSumKey(triggersumkey);
  uint16_t tmp = (triggersumkey >> kBitShiftPrimitiveLocId) & 0x1ffU;

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    if (detId == TriggerDefs::DetectorId::hcalinDId ||
        detId == TriggerDefs::DetectorId::hcaloutDId ||
        detId == TriggerDefs::DetectorId::hcalDId)
    {
      return tmp % 3;
    }
    if (detId == TriggerDefs::DetectorId::emcalDId)
    {
      return tmp % 12;
    }
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return 0;
  }

  return UINT16_MAX;
}

uint16_t
TriggerDefs::getSumLocId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint16_t tmp = (triggersumkey >> kBitShiftSumLocId) & 0x1fU;
  return tmp;
}
uint16_t
TriggerDefs::getSumPhiId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t primId = getPrimitiveId_from_TriggerSumKey(triggersumkey);
  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    return tmp / 4;
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return tmp / 12;
  }

  return UINT16_MAX;
}
uint16_t
TriggerDefs::getSumEtaId(const TriggerDefs::TriggerSumKey triggersumkey)
{
  uint32_t primId = getPrimitiveId_from_TriggerSumKey(triggersumkey);
  uint16_t tmp = TriggerDefs::getSumLocId(triggersumkey);

  if (primId == TriggerDefs::PrimitiveId::calPId)
  {
    return tmp % 4;
  }
  if (primId == TriggerDefs::PrimitiveId::jetPId)
  {
    return tmp % 12;
  }

  return UINT16_MAX;
}

uint32_t
TriggerDefs::GetTowerInfoKey(const TriggerDefs::DetectorId detId, const uint16_t iprim, const uint16_t isum, const uint16_t itower)
{
  unsigned int phibin = 0;
  unsigned int etabin = 0;
  if (detId == TriggerDefs::DetectorId::emcalDId)
  {
    phibin = 8 * (iprim / 12) + 2 * (isum / 4) + (itower / 2);
    etabin = 8 * (iprim % 12) + 2 * (isum % 4) + (itower % 2);
    return TowerInfoDefs::encode_emcal(etabin, phibin);
    ;
  }
  if (detId == TriggerDefs::DetectorId::hcalinDId ||
      detId == TriggerDefs::DetectorId::hcaloutDId ||
      detId == TriggerDefs::DetectorId::hcalDId)
  {
    etabin = (iprim % 3) * 8 + (isum % 4) * 2 + itower % 2;
    phibin = (iprim / 3) * 8 + (isum / 4) * 2 + itower / 2;
    return TowerInfoDefs::encode_hcal(etabin, phibin);
  }

  return UINT16_MAX;
}
TriggerDefs::TriggerId TriggerDefs::GetTriggerId(const std::string& trigger)
{
  if (trigger == "NONE")
  {
    return TriggerDefs::TriggerId::noneTId;
  }
  if (trigger == "MBD")
  {
    return TriggerDefs::TriggerId::mbdTId;
  }
  if (trigger == "JET")
  {
    return TriggerDefs::TriggerId::jetTId;
  }
  if (trigger == "PHOTON")
  {
    return TriggerDefs::TriggerId::photonTId;
  }
  if (trigger == "PAIR")
  {
    return TriggerDefs::TriggerId::pairTId;
  }
  if (trigger == "COSMIC")
  {
    return TriggerDefs::TriggerId::cosmicTId;
  }
  if (trigger == "COSMIC_COIN")
  {
    return TriggerDefs::TriggerId::cosmic_coinTId;
  }
  if (trigger == "PHYSICS")
  {
    return TriggerDefs::TriggerId::physicsTId;
  }

  return TriggerDefs::TriggerId::noneTId;
}

TriggerDefs::DetectorId TriggerDefs::GetDetectorId(const std::string& detector)
{
  if (detector == "NONE")
  {
    return TriggerDefs::DetectorId::noneDId;
  }
  if (detector == "MBD")
  {
    return TriggerDefs::DetectorId::mbdDId;
  }
  if (detector == "HCALIN")
  {
    return TriggerDefs::DetectorId::hcalinDId;
  }
  if (detector == "HCALOUT")
  {
    return TriggerDefs::DetectorId::hcaloutDId;
  }
  if (detector == "EMCAL")
  {
    return TriggerDefs::DetectorId::emcalDId;
  }
  if (detector == "CAL")
  {
    return TriggerDefs::DetectorId::calDId;
  }
  if (detector == "HCAL")
  {
    return TriggerDefs::DetectorId::hcalDId;
  }

  return TriggerDefs::DetectorId::noneDId;
}
TriggerDefs::PrimitiveId TriggerDefs::GetPrimitiveId(const std::string& primitive)
{
  if (primitive == "NONE")
  {
    return TriggerDefs::PrimitiveId::nonePId;
  }
  if (primitive == "MBD")
  {
    return TriggerDefs::PrimitiveId::mbdPId;
  }
  if (primitive == "HCALIN" ||
      primitive == "HCALOUT" ||
      primitive == "HCAL" ||
      primitive == "EMCAL")
  {
    return TriggerDefs::PrimitiveId::calPId;
  }
  if (primitive == "JET")
  {
    return TriggerDefs::PrimitiveId::jetPId;
  }
  if (primitive == "PAIR")
  {
    return TriggerDefs::PrimitiveId::pairPId;
  }

  return TriggerDefs::PrimitiveId::nonePId;
}
