#ifndef VARARRAY_VARIABLEARRAY_H
#define VARARRAY_VARIABLEARRAY_H

#include <phool/PHObject.h>

#include <vector>

class VariableArray : public PHObject
{
 public:
  VariableArray(const unsigned int idval = 0);
  virtual ~VariableArray();

  void identify(std::ostream &os=std::cout) const;

  // Here are the very explicit set routines...
  void set_val(const std::vector<short> &vec); 
  const short int *get_array() const {return sval;}
  unsigned int get_array_size() const {return nVal;}
  int Id() const {return id;}
  void Reset();

 protected:
  int id;
  unsigned int nVal;
  short *sval; //[nVal]

  ClassDef(VariableArray,1)
};

#endif /* VARIABLEARRAY */

