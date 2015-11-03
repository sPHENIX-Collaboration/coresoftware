//  Declaration of class PdbADCChan
//  Purpose: Example of a calibration object
//  Author: Matthias Messer

#ifndef __PDBADCCHAN_HH__
#define __PDBADCCHAN_HH__

#include "PdbCalChan.hh"

class PdbADCChan : public PdbCalChan 
{
public:
   PdbADCChan();
  virtual ~PdbADCChan( ){}

   size_t getNDim() const { return nDim; }
   const char* getParName(const size_t) const;

   float getLowGain()     const {return ADCParameter[0];}
   float getHighGain()    const {return ADCParameter[1];}
   float getConvert()     const {return ADCParameter[2];}
   float getLowGainErr()  const {return ADCParError[0];}
   float getHighGainErr() const {return ADCParError[1];}
   float getConvertErr()  const {return ADCParError[2];}
   
   float getParameter(const size_t) const;
   float getParError(const size_t) const;
   
   void setLowGain(const float val)     { ADCParameter[0] = val;}
   void setHighGain(const float val)    { ADCParameter[1] = val;}
   void setConvert(const float val)     { ADCParameter[2] = val;}
   void setLowGainErr(const float val)  { ADCParError[0] = val;}
   void setHighGainErr(const float val) { ADCParError[1] = val;}
   void setConvertErr(const float val)  { ADCParError[2] = val;}

   void setParameter(const size_t, const float);
   void setParError(const size_t, const float);
   
   virtual void print() const;

private:
   void zero();
   
private:
   size_t nDim;
   float ADCParameter[3];
   float ADCParError[3];

  ClassDef(PdbADCChan,1);
};

#endif /* __PDBADCCHAN_HH__ */
