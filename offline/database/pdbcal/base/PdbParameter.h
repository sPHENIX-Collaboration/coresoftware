//  Declaration of class PdbParameter
//  Purpose: single parameter storage class 
//  Author: federica

#ifndef __PDBPARAMETER_HH__
#define __PDBPARAMETER_HH__

#include "PdbCalChan.h"

class PdbParameter : public PdbCalChan {
public:
  PdbParameter();
  PdbParameter(const float);
  PdbParameter(const float, const char* name); 
  virtual ~PdbParameter() {}

  float getParameter() const  { return thePar;  }
  const char* getName() const { return theName; }

  void  setParameter(const float val) { thePar = val; }
  void  setName(const char* name);

  virtual void print() const;

private:

  float thePar;
  char  theName[20];

  ClassDef(PdbParameter,1);
};

#endif /* __PDBPARAMETER_DDL__ */
