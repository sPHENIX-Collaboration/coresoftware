//  Declaration of class PdbParameter
//  Purpose: single parameter storage class 
//  Author: federica

#ifndef PDBCAL_BASE_PDBPARAMETER_H
#define PDBCAL_BASE_PDBPARAMETER_H

#include "PdbCalChan.h"

#include <string>

class PdbParameter : public PdbCalChan {
public:
  PdbParameter(); // this ctor should not be called but it cannot be 
                  // made private since CINT needs a 
                  // default ctor when reading from file

  PdbParameter(const double, const std::string &name); 
  ~PdbParameter() override {}

  double getParameter() const  { return thePar;  }
  const std::string getName() const { return theName; }

  void  setParameter(const double val) { thePar = val; }
  void  setName(const std::string &name) {theName = name;}

  void print() const override;

protected:

  double thePar;
  std::string  theName;

  ClassDefOverride(PdbParameter,1);
};

#endif /* PDBCAL_BASE_PDBPARAMETER_H */
