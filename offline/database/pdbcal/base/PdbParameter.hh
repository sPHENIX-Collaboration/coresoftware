//  Declaration of class PdbParameter
//  Purpose: single parameter storage class 
//  Author: federica

#ifndef __PDBPARAMETER_HH__
#define __PDBPARAMETER_HH__

#include "PdbCalChan.hh"

#include <string>

class PdbParameter : public PdbCalChan {
public:
  PdbParameter();
  PdbParameter(const float);
  PdbParameter(const float, const char* name); 
  virtual ~PdbParameter() {}

  float getParameter() const  { return thePar;  }
  const std::string getName() const { return theName; }

  void  setParameter(const float val) { thePar = val; }
  void  setName(const std::string &name) {theName = name;}

  virtual void print() const;

private:

  double thePar;
  std::string  theName;

  ClassDef(PdbParameter,1);
};

#endif /* __PDBPARAMETER_DDL__ */
