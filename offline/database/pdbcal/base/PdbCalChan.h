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
#ifndef PDBCAL_BASE_PDBCALCHAN_H
#define PDBCAL_BASE_PDBCALCHAN_H

#include <phool/PHObject.h>

class PdbCalChan : public PHObject
{
 public:
  PdbCalChan() {}
  virtual ~PdbCalChan() {}

  virtual void print() const = 0;

  ClassDef(PdbCalChan, 1);
};

#endif /* PDBCAL_BASE_PDBCALCHAN_H */
