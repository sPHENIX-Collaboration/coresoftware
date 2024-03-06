
#ifndef __TRIGGERPRIMITIVECONTAINERV1_H
#define __TRIGGERPRIMITIVECONTAINERV1_H

#include <string>
#include <ostream>
#include <iostream>
#include <phool/PHObject.h>
#include "TriggerPrimitiveContainer.h"
#include <map>
#include <utility>

//#include <TClonesArray.h>
//#include <calobase/RawTowerDefs.h> 
///

class TriggerPrimitiveContainerv1 : public TriggerPrimitiveContainer
{
 public:
  ///
  TriggerPrimitiveContainerv1();
  ///
  ~TriggerPrimitiveContainerv1();

  /// Clear Event from memory
  virtual void Reset() override;
  virtual void identify(std::ostream& os = std::cout) const override;
  virtual int isValid() const override;

  virtual void setTriggerType(TriggerDefs::TriggerId triggerid) {m_triggerkey = TriggerDefs::getTriggerKey(triggerid);}

  virtual TriggerPrimitive* get_primitive_at_key(TriggerDefs::TriggerPrimKey /* index */ ) override;

  virtual void add_primitive(TriggerDefs::TriggerPrimKey , TriggerPrimitive* ) override;

  virtual TriggerDefs::TriggerKey getTriggerKey() {return m_triggerkey;}

  ConstRange getTriggerPrimitives() const;  
  Range getTriggerPrimitives();

  virtual size_t size() override { return _primitives.size();}

 protected:

  TriggerDefs::TriggerKey m_triggerkey = TriggerDefs::TRIGGERKEYMAX;
  Map _primitives;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerPrimitiveContainerv1,1);
};

#endif
