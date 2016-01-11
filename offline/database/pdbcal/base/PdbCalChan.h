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
#ifndef PDBCALCHAN_HH__
#define PDBCALCHAN_HH__

#include <phool/PHObject.h>

class PdbCalChan : public PHObject{

public:
  PdbCalChan() {}
  virtual ~PdbCalChan() {}
  
  virtual void print() const = 0;

  ClassDef(PdbCalChan,1);
};

#endif /* PDBCALCHAN_HH__ */
