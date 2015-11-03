#ifndef __PDBFLOATVECTOR_HH__
#define __PDBFLOATVECTOR_HH__

#include "PdbCalChan.hh"
#include <vector>

class PdbFloatVector : public PdbCalChan
{
  public:
  PdbFloatVector(){}
  virtual ~PdbFloatVector();

  void add_float(float fl);
  virtual void print() const;
  
  std::vector<float> getVector() const  { return _floatvec; }
  float getValue(int pos) const;
  int getLength() const {return _floatvec.size(); }

private:
  
  std::vector<float> _floatvec;
  
  ClassDef(PdbFloatVector,1);
};

#endif /* __PDBFLOATVECTOR_HH__ */
