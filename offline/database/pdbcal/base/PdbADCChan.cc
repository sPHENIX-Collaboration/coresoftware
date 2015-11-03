//  Implementation of class PdbADCChan
//  Author: Matthias Messer

#include <iostream>
#include "PdbADCChan.hh"

using namespace std;

PdbADCChan::PdbADCChan() : nDim(3)
{
   zero();
}

void 
PdbADCChan::zero()
{
  for (size_t i = 0; i < nDim; i++) 
    {
      ADCParameter[i] = 0.0;
      ADCParError[i]  = 0.0;
    }
}

float 
PdbADCChan::getParameter(const size_t i) const
{
  // Returns the parameter value at a given index location, 0.0 if the index is out of range.

  if (i < nDim)
    {
      return ADCParameter[i];
    }
  
  return 0.0;
}

float 
PdbADCChan::getParError(const size_t i) const
{
   // Returns the parameter error at a given index location, 0.0 if the index is out of range.

  if (i < nDim)
    {
      return ADCParError[i];
    }

  return 0.0;
}

const char* 
PdbADCChan::getParName(const size_t i) const
{
   // Returns the parameter error at a given index location. Returns 0.0 if the index is out of range.

  switch (i) 
    {
    case 0:
      return "Low Gain";
    case 1:
      return "High Gain";
    case 2:
      return "Conversion Factor";
    default:
      return 0;
    }
}

void 
PdbADCChan::setParameter(const size_t i, const float val)
{
  if (i < nDim)
    {
      ADCParameter[i] = val;
    }
   else
     {
      cout << "PdbADCChan::SetParameter - Index = " << i 
	   << " is out of range. [0.." << nDim-1 << "] is valid." 
	   << endl;
     }
}

void 
PdbADCChan::setParError(const size_t i, const float val)
{
   if (i < nDim)
     {
       ADCParError[i] = val;
     }
   else
     {
       cout << "PdbADCChan::SetParError - Index = " << i  
	    << " is out of range. [0.." << nDim-1 << "] is valid." 
	    << endl;
     }
}

void 
PdbADCChan::print() const
{
   cout << "LowGain  = " << ADCParameter[0] << " +- " << ADCParError[0] << endl
	<< "HighGain = " << ADCParameter[1] << " +- " << ADCParError[1] << endl
	<< "Convert  = " << ADCParameter[2] << " +- " << ADCParError[2] << endl;
}
