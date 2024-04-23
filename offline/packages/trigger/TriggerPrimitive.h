#ifndef TRIGGER_TRIGGERPRIMITIVE_H
#define TRIGGER_TRIGGERPRIMITIVE_H

#include "TriggerDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <ostream>
#include <vector>

///
class TriggerPrimitive : public PHObject
{
 public:
  typedef std::map<TriggerDefs::TriggerSumKey, std::vector<unsigned int>*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  TriggerPrimitive() = default;
  virtual ~TriggerPrimitive() override = default;

  /// Clear Event from memory
  virtual void Reset() override { return; }
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override { return 0; }

  virtual std::vector<unsigned int>* get_sum_at_key(TriggerDefs::TriggerSumKey) { return nullptr; }

  virtual void add_sum(TriggerDefs::TriggerSumKey /*key*/, std::vector<unsigned int>* /*sum*/) { return; }

  virtual ConstRange getSums() const;
  virtual Range getSums();

  virtual size_t size() { return 0; }

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitive, 1);
};

#endif
