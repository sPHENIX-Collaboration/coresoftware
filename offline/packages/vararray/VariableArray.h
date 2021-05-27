#ifndef VARARRAY_VARIABLEARRAY_H
#define VARARRAY_VARIABLEARRAY_H

#include <phool/PHObject.h>

#include <iostream>
#include <vector>

class VariableArray : public PHObject
{
 public:
  VariableArray(const unsigned int idval = 0);
  ~VariableArray() override;

  void identify(std::ostream &os = std::cout) const override;

  // Here are the very explicit set routines...
  void set_val(const std::vector<short> &vec);
  const short int *get_array() const { return sval; }
  unsigned int get_array_size() const { return nVal; }
  int Id() const { return id; }
  void Reset() override;

 protected:
  int id;
  unsigned int nVal;
  short *sval;  //[nVal]

  ClassDefOverride(VariableArray, 1)
};

#endif /* VARIABLEARRAY */
