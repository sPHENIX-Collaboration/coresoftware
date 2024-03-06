#ifndef __LL1OUT_H
#define __LL1OUT_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <ostream>
#include <phool/PHObject.h>


///
class LL1Out : public PHObject
{
 public:
  
  ///
  LL1Out();
  ///
  virtual ~LL1Out();
  
  /// Clear Event from memory
  virtual void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream 
  */
  virtual void identify(std::ostream& os = std::cout);

  /// isValid returns non zero if object contains vailid data
  virtual int isValid();

  // Get Trigger Type
  //  virtual std::string get_TriggerType() const {return _trigger_type;}

  // Set Trigger Type  
  //virtual void set_TriggerType(std::string &triggertype) const { _trigger_type = triggertype;}

 protected:
  
  virtual void Init();
  
 private:
  
  std::string _trigger_type;

 private: // so the ClassDef does not show up with doc++
  ClassDefOverride(LL1Out,1);
};

#endif
