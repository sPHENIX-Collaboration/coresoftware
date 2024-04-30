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

  TriggerPrimitiveContainerv1(const std::string& triggertype);
  ///
  ~TriggerPrimitiveContainerv1() override;

  /// Clear Event from memory
  void Reset() override;
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;

  void setTriggerType(TriggerDefs::TriggerId triggerid) override { m_triggerkey = TriggerDefs::getTriggerKey(triggerid); }

  TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey /* index */) override;

  void add_primitive(TriggerDefs::TriggerPrimKey, TriggerPrimitive*) override;

  TriggerDefs::TriggerKey getTriggerKey() { return m_triggerkey; }

  ConstRange getTriggerPrimitives() const override;
  Range getTriggerPrimitives() override;

  size_t size() override { return _primitives.size(); }

 private:  // so the ClassDef does not show up with doc++
  TriggerDefs::TriggerKey m_triggerkey = TriggerDefs::TRIGGERKEYMAX;
  Map _primitives;

  ClassDefOverride(TriggerPrimitiveContainerv1, 1);
};

#endif
