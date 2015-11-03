#ifndef __PDBDOUBLE_HH__
#define __PDBDOUBLE_HH__

#include "PdbCalChan.hh"

class PdbDouble : public PdbCalChan 
{
public:
   PdbDouble();
   ~PdbDouble();

   double getValue() const {return TheValue;}
   void setValue(double x) {TheValue=x;}

   virtual void print() const;

private:
  double TheValue;

  ClassDef(PdbDouble,1);

};

#endif /* __PDBDOUBLE_HH__ */
