//  Declaration of class PdbParameterError
//  Purpose: single parameter storage class
//  Author: Cesar & federica

#ifndef PDBCAL_BASE_PDBPARAMETERERROR_H
#define PDBCAL_BASE_PDBPARAMETERERROR_H

#include "PdbParameter.h"

class PdbParameterError : public PdbParameter
{
 public:
  PdbParameterError();
  PdbParameterError(const double, const double, const std::string &name);
  virtual ~PdbParameterError() {}

  double getParameterError() const { return theParError; }

  void setParameterError(const double val) { theParError = val; }

  virtual void print() const;

 protected:
  double theParError;

  ClassDef(PdbParameterError, 1);
};

#endif /* PDBCAL_BASE_PDBPARAMETERERROR_H */
