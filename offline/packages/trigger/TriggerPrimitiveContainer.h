#ifndef __TRIGGERPRIMITIVECONTAINER_H
#define __TRIGGERPRIMITIVECONTAINER_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <ostream>
#include <phool/PHObject.h>
#include <map>
#include <climits>
#include "TriggerPrimitive.h"
#include "TriggerDefs.h"

///
class TriggerPrimitiveContainer : public PHObject
{
 public:
  TriggerPrimitiveContainer();
  virtual ~TriggerPrimitiveContainer();
  
  /// Clear Event from memory
  virtual void Reset() override;
  void identify(std::ostream& os = std::cout) const override;
  int isValid() const override;

  virtual TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey) { return nullptr; }

  virtual void add_primitive(TriggerDefs::TriggerKey key, TriggerPrimitive* primitive);

  virtual size_t size() {return 0;}

 private:
  
  TriggerDefs::TriggerKey m_triggerkey = TriggerDefs::TRIGGERKEYMAX;
  //  Map _primitives;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitiveContainer,1);
};

#endif
