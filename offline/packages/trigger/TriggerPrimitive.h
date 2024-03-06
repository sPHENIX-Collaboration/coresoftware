#ifndef __TRIGGERPRIMITIVE_H__
#define __TRIGGERPRIMITIVE_H__

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <ostream>
#include <phool/PHObject.h>
#include <map>
#include <climits>
#include <vector>
#include "TriggerDefs.h"

///
class TriggerPrimitive : public PHObject
{
 public:
  typedef std::map<TriggerDefs::TriggerSumKey, std::vector<unsigned int>*> Map;
  typedef Map::const_iterator ConstIter;
  typedef Map::iterator Iter;
  typedef std::pair<Iter, Iter> Range;
  typedef std::pair<ConstIter, ConstIter> ConstRange;

  TriggerPrimitive();
  TriggerPrimitive(TriggerDefs::TriggerPrimKey key);
  virtual ~TriggerPrimitive();
  
  /// Clear Event from memory
  virtual void Reset() override;
  virtual void identify(std::ostream& os = std::cout);
  virtual int isValid();

  virtual std::vector<unsigned int>* get_sum_at_key(TriggerDefs::TriggerSumKey );

  virtual void add_sum(TriggerDefs::TriggerSumKey key, std::vector<unsigned int> *sum);

  virtual TriggerDefs::TriggerPrimKey getTriggerPrimitiveKey() { return m_triggerprimkey;}

  virtual ConstRange getSums() const;
  virtual Range getSums();
  
  virtual size_t size() {return _sums.size();}

 private:
  
  TriggerDefs::TriggerPrimKey m_triggerprimkey = TriggerDefs::TRIGGERPRIMKEYMAX;
  Map _sums;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitive,1);
};

#endif
