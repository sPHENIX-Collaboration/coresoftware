//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbCalChan
//
//  Purpose: Calibration channel base class
//
//  Description:
//
//  Author: Matthias Messer
//-----------------------------------------------------------------------------
#ifndef __PDBCALCHAN_HH__
#define __PDBCALCHAN_HH__

#include <TObject.h>

class PdbCalChan : public TObject{

public:
  PdbCalChan() {}
  virtual ~PdbCalChan() {}
  
  virtual void print() const = 0;

  ClassDef(PdbCalChan,1);
};

#endif /* __PDBCALCHAN_HH__ */
