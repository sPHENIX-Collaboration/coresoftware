//  Declaration of class PdbIndex
//  Purpose: User defined storage class
//  Author: federica

#ifndef __PDBINDEX_DDL__
#define __PDBINDEX_DDL__

#include "PdbCalChan.h"

#include <phool/phool.h>

class PdbIndex : public PdbCalChan 
{
public:

  PdbIndex();
  PdbIndex(const int mini, const int maxi);
  PdbIndex(const int mini, const int maxi, const char* name);
  PdbIndex(const int mini, const int maxi, const int val, const char* name); 
  
  virtual ~PdbIndex(){}

  int operator()() {return theValue;}

  virtual void print() const;
  
  int   getMin()      const { return theMinimum;}
  int   getMax()      const { return theMaximum;}
  int   getValue()    const { return theValue;}
  int   getNumberOf() const { return (theMaximum - theMinimum +1);}
  
  void setMin(const int val)        { theMinimum = val; }
  void setMax(const int val)        { theMaximum = val; }
  void setName(const char* name);
  PHBoolean setValue(const int val);

private:
  int   theMinimum;
  int   theMaximum;
  int   theValue;
  char  theName[20];

  ClassDef(PdbIndex,1);
};

#endif /* __PDBINDEX_DDL__ */
