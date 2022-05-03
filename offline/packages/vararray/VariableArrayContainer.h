#ifndef VARARRAY_VARIABLEARRAYCONTAINER_H
#define VARARRAY_VARIABLEARRAYCONTAINER_H

#include <phool/PHObject.h>

#include <iostream>

class TObjArray;
class VariableArray;

class VariableArrayContainer : public PHObject
{
 public:
  VariableArrayContainer();
  ~VariableArrayContainer() override;

  void identify(std::ostream &os = std::cout) const override;
  void AddVarArray(VariableArray *var);
  // Here are the very explicit set routines...
  void Reset() override;

 protected:
  TObjArray *arraycontainer;

  ClassDefOverride(VariableArrayContainer, 1)
};

#endif
