//  Implementation of class PdbCoordinate
//  Author: messer

#include <iostream>
#include "PdbCoordinate.hh"

using namespace std;

PdbCoordinate::PdbCoordinate()
{
  nDim = 3;
  for (size_t ix = 0; ix < nDim; ix++) 
    {
      fCoordinate[ix] = 0.0;
      fCoordinateError[ix] = 0.0;
    }
}

PdbCoordinate::PdbCoordinate(const PdbCoordinate &rhs) 
{
  nDim = 3;
  for (size_t ix = 0; ix < nDim; ix++) {
    fCoordinate[ix] = rhs.fCoordinate[ix];
    fCoordinateError[ix] = rhs.fCoordinateError[ix];
  }

}

PdbCoordinate::PdbCoordinate(float xx, float yy, float zz)
{
  nDim = 3;
  fCoordinate[0] = xx;
  fCoordinate[1] = yy;
  fCoordinate[2] = zz;
  for (size_t ix = 0; ix < nDim; ix++) {
    fCoordinateError[ix] = 0.0;
  }
}

PdbCoordinate::~PdbCoordinate()
{
}

float PdbCoordinate::getParameter(size_t index) const
{
   //
   // Returns the parameter value at a given index location. Returns 0.0 if the index is out of range.
   //
   switch(index) {
   case x:
      return fCoordinate[x];
   case y:
      return fCoordinate[y];
   case z:
      return fCoordinate[z];
   default:
      return 0.0;
   }
}



float PdbCoordinate::getParError(size_t index) const
{
   //
   // Returns the parameter error at a given index location. Returns 0.0 if the index is out of range.
   //
   switch(index) {
   case x:
      return fCoordinateError[x];
   case y:
      return fCoordinateError[y];
   case z:
      return fCoordinateError[z];
   default:
      return 0.0;
   }
}

const char* PdbCoordinate::getParName(size_t index) const
{
   //
   // Returns the parameter error at a given index location. Returns 0.0 if the index is out of range.
   //
   switch(index) {
   case x:
      return "X";
   case y:
      return "Y";
   case z:
      return "Z";
   default:
      return 0;
   }
}

void PdbCoordinate::setParameter(size_t index,float theParValue)
{
  switch (index) {
  case x:
    fCoordinate[x]=theParValue;
    break;
  case y:
    fCoordinate[y]=theParValue;
    break;
  case z:
    fCoordinate[z]=theParValue;
    break;
  default:
    cout << "PdbCoordinate::setParameter - Index value = " 
	 << index  << " is out of range. [0.."
	 << nDim-1 << "] is valid." << endl;
  }
}

void PdbCoordinate::setAllParameters(float first, float second, float third)
{
    fCoordinate[x]=first;
    fCoordinate[y]=second;
    fCoordinate[z]=third;
}



void PdbCoordinate::setParError(size_t index, float theParError)
{
  switch (index) {
  case x:
    fCoordinateError[x]=theParError;
    break;
  case y:
    fCoordinateError[y]=theParError;
    break;
  case z:
    fCoordinateError[z]=theParError;
    break;
  default:
    cout << "PdbCoordinate::setParError - Index value = " 
	 << index  << " is out of range. [0.."
	 << nDim-1 << "] is valid." << endl;
  }
}

void PdbCoordinate::setAllParErrors(float first, float second, float third)
{
   fCoordinateError[x]=first;
   fCoordinateError[y]=second;
   fCoordinateError[z]=third; 
}

void 
PdbCoordinate::print() const
{ 
  for (size_t i =0; i<nDim; i++)
    {
      cout << getParName(i) << "\t\t" << fCoordinate[i] << "\t+/-\t" << fCoordinateError[i] << endl;
    }
}

