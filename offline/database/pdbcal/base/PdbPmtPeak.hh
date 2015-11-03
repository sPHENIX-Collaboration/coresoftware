//-----------------------------------------------------------------------------
//
//  The pdbcal package
//  Copyright (C) PHENIX collaboration, 1999
//
//  Declaration of class PdbPmtPeak
//
//  Purpose: User defined storage class
//
//  Description:
//
//  Author: ohnishi
//-----------------------------------------------------------------------------
#ifndef __PDBPMTPEAK_HH__
#define __PDBPMTPEAK_HH__

#include "PdbCalChan.hh"

class PdbPmtPeak : public PdbCalChan {
 public:
  PdbPmtPeak();
  virtual ~PdbPmtPeak(){}
  
  void setPeakChannel(const float peak ) {PeakChannel=peak;}
  void setDeviation(const float devi ) {Deviation=devi;} 
  void setStatus(const int stat ) {Status=stat;}
  
  float getPeakChannel() const {return PeakChannel;}
  float getDeviation()const  {return Deviation;} 
  int   getStatus() const {return Status;}
  
  virtual void print() const;
  
private:
  float PeakChannel;
  float Deviation;
  int   Status;

  ClassDef(PdbPmtPeak,1);

};

#endif /* __PDBPMTPEAK_HH__ */
