#ifndef TRIGGER_TRIGGERPRIMITIVECONTAINERV1_H
#define TRIGGER_TRIGGERPRIMITIVECONTAINERV1_H

#include "TriggerPrimitive.h"
#include "TriggerPrimitiveContainer.h"
#include "TriggerPrimitivev1.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <utility>

///

class TriggerPrimitiveContainerv1 : public TriggerPrimitiveContainer
{
 public:
  typedef std::map<TriggerDefs::TriggerPrimKey, TriggerPrimitive*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  ///
  TriggerPrimitiveContainerv1() = default;

  TriggerPrimitiveContainerv1(TriggerDefs::TriggerId tid, TriggerDefs::DetectorId did);
  ///
  ~TriggerPrimitiveContainerv1() override;

  /// Clear Event from memory
  void Reset() override;
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;

  void setTriggerId(TriggerDefs::TriggerId triggerid) override { m_triggerid = triggerid ; };
  void setDetectorId(TriggerDefs::DetectorId detectorid) override { m_detectorid = detectorid ; };
  void setPrimitiveId(TriggerDefs::PrimitiveId primitiveid) override { m_primitiveid = primitiveid ; };

  TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey /* index */) override;

  void add_primitive(TriggerDefs::TriggerPrimKey, TriggerPrimitive*) override;
  TriggerDefs::TriggerId getTriggerId() { return m_triggerid; }
  TriggerDefs::DetectorId getDetectorId() { return m_detectorid; }
  TriggerDefs::PrimitiveId getPrimitiveId() { return m_primitiveid; }

  ConstRange getTriggerPrimitives() const override;
  Range getTriggerPrimitives() override;

  size_t size() override { return _primitives.size(); }

 private:  // so the ClassDef does not show up with doc++
  TriggerDefs::TriggerId m_triggerid = TriggerDefs::TriggerId::noneTId;
  TriggerDefs::DetectorId m_detectorid = TriggerDefs::DetectorId::noneDId;
  TriggerDefs::PrimitiveId m_primitiveid = TriggerDefs::PrimitiveId::nonePId;

  Map _primitives;

  ClassDefOverride(TriggerPrimitiveContainerv1, 1);
};

#endif
