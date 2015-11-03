//-----------------------------------------------------------------------------
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbPmtFitPar
//
//  Purpose: User defined storage class
//
//  Description:
//
//  Author: ohnishi
//-----------------------------------------------------------------------------
#ifndef __PDBPMTFITPAR_HH__
#define __PDBPMTFITPAR_HH__

#include "PdbCalChan.hh"

class PdbPmtFitPar : public PdbCalChan {
public:
	  PdbPmtFitPar();
	  virtual ~PdbPmtFitPar();

          void setPar0(float par){Par0=par;}
          void setPar1(float par){Par1=par;}
          void setPar2(float par){Par2=par;}
          void setPar3(float par){Par3=par;}
          void setPar4(float par){Par4=par;}
          void setChi2(float c2) {Chi2=c2;}
          void setStatus( int stat ) {Status=stat;}

          float getPar0() const {return Par0;}
          float getPar1() const {return Par1;}
          float getPar2() const {return Par2;}
          float getPar3() const {return Par3;}
          float getPar4() const {return Par4;}
          float getChi2() const {return Chi2;}
          int   getStatus() const {return Status;}

  virtual void print() const;
  PdbPmtFitPar& operator = (const PdbPmtFitPar &p);

private:
          float Par0;
          float Par1;
          float Par2;
          float Par3;
          float Par4;
          float Chi2;
          int   Status;

  ClassDef(PdbPmtFitPar,1);

};

#endif /* __PDBPMTFITPAR_HH__ */
