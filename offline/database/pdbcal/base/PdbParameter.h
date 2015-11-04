//  Declaration of class PdbParameter
//  Purpose: single parameter storage class 
//  Author: federica

#ifndef PDBPARAMETER_HH__
#define PDBPARAMETER_HH__

#include "PdbCalChan.h"

#include <string>

class PdbParameter : public PdbCalChan {
public:
  PdbParameter(); // this ctor should not be called but it cannot be 
                  // made private since CINT needs a 
                  // default ctor when reading from file

  PdbParameter(const double, const std::string &name); 
  virtual ~PdbParameter() {}

  double getParameter() const  { return thePar;  }
  const std::string getName() const { return theName; }

  void  setParameter(const double val) { thePar = val; }
  void  setName(const std::string &name) {theName = name;}

  virtual void print() const;

protected:

  double thePar;
  std::string  theName;

  ClassDef(PdbParameter,1);
};

#endif /* __PDBPARAMETER_DDL__ */
