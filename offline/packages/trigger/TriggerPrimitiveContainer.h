#ifndef TRIGGER_TRIGGERPRIMITIVECONTAINER_H
#define TRIGGER_TRIGGERPRIMITIVECONTAINER_H

#include "TriggerDefs.h"
#include "TriggerPrimitive.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <ostream>

///
class TriggerPrimitiveContainer : public PHObject
{
 public:
  typedef std::map<TriggerDefs::TriggerPrimKey, TriggerPrimitive*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  TriggerPrimitiveContainer() = default;
  virtual ~TriggerPrimitiveContainer() override = default;

  /// Clear Event from memory
  virtual void Reset() override { return; };
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 1; }
  virtual void setTriggerType(TriggerDefs::TriggerId /*triggerid*/) { return; }
  virtual TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey) { return nullptr; }

  virtual void add_primitive(TriggerDefs::TriggerPrimKey, TriggerPrimitive*) { return; }

  virtual size_t size() { return 0; }

  virtual ConstRange getTriggerPrimitives() const;
  virtual Range getTriggerPrimitives();

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitiveContainer, 1);
};

#endif
