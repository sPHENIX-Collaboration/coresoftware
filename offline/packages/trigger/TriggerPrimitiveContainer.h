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
  typedef std::map<TriggerDefs::TriggerPrimKey, TriggerPrimitive*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;


  TriggerPrimitiveContainer();
  virtual ~TriggerPrimitiveContainer();
  
  /// Clear Event from memory
  virtual void Reset() override;
  virtual void identify(std::ostream& os = std::cout);
  virtual int isValid();

  virtual TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey) { return nullptr; }

  virtual void add_primitive(TriggerDefs::TriggerKey key, TriggerPrimitive* primitive);

  virtual size_t size() {return 0;}

 private:
  
  TriggerDefs::TriggerKey m_triggerkey = TriggerDefs::TRIGGERKEYMAX;
  Map _primitives;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitiveContainer,1);
};

#endif
