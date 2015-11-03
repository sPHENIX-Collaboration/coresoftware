//  Declaration of class PdbCoordinate
//  Purpose: User defined storage class
//  Author: messer

#ifndef __PDBCOORDINATE_DDL__
#define __PDBCOORDINATE_DDL__

#include "PdbCalChan.hh"
#ifndef __CINT__
#include <cstddef>
#endif

class PdbCoordinate : public PdbCalChan 
{
public:
   enum {x=0, y=1, z=2};
   
   PdbCoordinate();
   PdbCoordinate(float x, float y, float z);
   PdbCoordinate(const PdbCoordinate &);
   virtual ~PdbCoordinate();

   friend PdbCoordinate operator- (const PdbCoordinate &,
				   const PdbCoordinate &);
   friend PdbCoordinate operator+ (const PdbCoordinate &,
				   const PdbCoordinate &);
   
   size_t getNdim() const { return nDim;}
   float  getParameter(size_t) const;
   float  getParError(size_t) const;
   const  char* getParName(size_t) const;

   void setParameter(size_t, float);
   void setParError(size_t, float);
   void setAllParameters(float, float, float);
   void setAllParErrors(float, float, float);

   virtual void print() const;

private:
   size_t nDim;
   float  fCoordinate[3];  
   float  fCoordinateError[3];   // Errors on the parameters.

  ClassDef(PdbCoordinate,1);

};

inline 
PdbCoordinate 
operator- (const PdbCoordinate & a, const PdbCoordinate & b) 
{
  size_t i;
  PdbCoordinate c;

  for (i = 0 ; i < c.getNdim(); i++) 
    {
      c.setParameter(i, a.getParameter(i) - b.getParameter(i));
    }

  return c;
}

inline 
PdbCoordinate 
operator+ (const PdbCoordinate & a, const PdbCoordinate & b) 
{
  size_t i;
  PdbCoordinate c;

  for (i = 0 ; i < c.getNdim(); i++) 
    {
      c.setParameter(i, a.getParameter(i) + b.getParameter(i));
    }
  
  return c;
}

#endif /* __PDBCOORDINATE_DDL__ */
