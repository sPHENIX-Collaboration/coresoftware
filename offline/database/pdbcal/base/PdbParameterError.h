//  Declaration of class PdbParameterError
//  Purpose: single parameter storage class 
//  Author: Cesar & federica

#ifndef PDBPARAMETERERROR_HH__
#define PDBPARAMETERERROR_HH__

#include "PdbParameter.h"

class PdbParameterError : public PdbParameter 
{
 public:
  PdbParameterError();
  PdbParameterError(const double, const double, const std::string &name); 
  virtual ~PdbParameterError() {}

  double getParameterError() const { return theParError; }

  void  setParameterError(const double val) { theParError = val; }

  virtual void print() const;

 protected:

  double theParError;

  ClassDef(PdbParameterError,1);
};

#endif /* PDBPARAMETERERROR_HH__ */
