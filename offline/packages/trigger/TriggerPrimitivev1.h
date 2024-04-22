#ifndef TRIGGER_TRIGGERPRIMITIVEV1_H__
#define TRIGGER_TRIGGERPRIMITIVEV1_H__

#include "TriggerDefs.h"
#include "TriggerPrimitive.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>
#include <ostream>
#include <vector>

///
class TriggerPrimitivev1 : public TriggerPrimitive
{
 public:
  TriggerPrimitivev1() = default;
  TriggerPrimitivev1(TriggerDefs::TriggerPrimKey key);
  virtual ~TriggerPrimitivev1() override = default;

  /// Clear Event from memory
  void Reset() override;
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;

  std::vector<unsigned int>* get_sum_at_key(TriggerDefs::TriggerSumKey) override;

  void add_sum(TriggerDefs::TriggerSumKey key, std::vector<unsigned int>* sum) override;

  ConstRange getSums() const override;
  Range getSums() override;

  size_t size() override { return _sums.size(); }

 private:
  TriggerDefs::TriggerPrimKey m_triggerprimkey = TriggerDefs::TRIGGERPRIMKEYMAX;
  Map _sums;

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitivev1, 1);
};

#endif
