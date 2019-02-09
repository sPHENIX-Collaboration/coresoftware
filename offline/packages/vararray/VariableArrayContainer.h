#ifndef VARARRAY_VARIABLEARRAYCONTAINER_H
#define VARARRAY_VARIABLEARRAYCONTAINER_H

#include <phool/PHObject.h>

#include <vector>

class TOjArray;
class VariableArray;

class VariableArrayContainer : public PHObject
{
 public:
  VariableArrayContainer();
  virtual ~VariableArrayContainer();

  void identify(std::ostream &os = std::cout) const;
  void AddVarArray(VariableArray *var);
  // Here are the very explicit set routines...
  void Reset();

 protected:
  TObjArray *arraycontainer;

  ClassDef(VariableArrayContainer, 1)
};

#endif
