//  Declaration of class PdbParameterError
//  Purpose: single parameter storage class 
//  Author: Cesar & federica

#ifndef __PDBPARAMETERERROR_DDL__
#define __PDBPARAMETERERROR_DDL__

#include <PdbCalChan.hh>

class PdbParameterError : public PdbCalChan {
public:
  PdbParameterError();
  PdbParameterError(const float);
  PdbParameterError(const float, const float);
  PdbParameterError(const float, const float, const char* name); 
  virtual ~PdbParameterError() {}

  float getParameter() const { return thePar; }
  float getParameterError() const { return theParError; }
  const char* getName() const { return theName; }

  void  setParameter(const float val) { thePar = val; }
  void  setParameterError(const float val) { theParError = val; }
  void  setName(const char* name);

  virtual void print() const;

private:

  float thePar;
  float theParError;
  char  theName[30];

  ClassDef(PdbParameterError,1);
};

#endif /* __PDBPARAMETERERROR_DDL__ */
