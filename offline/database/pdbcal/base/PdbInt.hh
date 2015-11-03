#ifndef __PDBINT_HH__
#define __PDBINT_HH__

#include "PdbCalChan.hh"

class PdbInt : public PdbCalChan 
{
public:
   PdbInt();
  virtual ~PdbInt() {}

   int getValue() const {return TheValue;}
   void setValue(const int x) {TheValue=x;}

   virtual void print() const;

private:
  int TheValue;

  ClassDef(PdbInt,1);

};

#endif /* __PDBINT_HH__ */
