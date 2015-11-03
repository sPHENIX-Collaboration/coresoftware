#include <iostream>
#include "PdbDouble.hh"

using namespace std;

PdbDouble::PdbDouble() 
{
  TheValue = 0;
}

PdbDouble::~PdbDouble()
{
}


void 
PdbDouble::print() const
{
  cout << TheValue << endl; 

}







