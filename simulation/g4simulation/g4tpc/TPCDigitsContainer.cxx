// TPC DIGITSCONTAINER class
// Stores TPC DIGITS per module
// notice that it stores a digits in vector and not in map
// no need of index access since reconstruction will expand per module
// i.e. stores compressed info
// features avoidance of ctor dtor calls
// Author: Carlos Perez

#include <map>
#include "TPCDigitsContainer.h"
#include "TPCDataTypes.h"
using namespace TPCDataTypes;

ClassImp(TPCDigitsContainer)

//=====
TPCDigitsContainer::TPCDigitsContainer()
{
  for(int i=0; i!=72; ++i)
    fNDigits[i] = 0;
}
//=====
TPCDigitsContainer::~TPCDigitsContainer()
{
  // freeing objects
  for(int i=0; i!=72; ++i)
    for(unsigned int j=0; j!=fDigits[i].size(); ++j)
      delete fDigits[i].at(j);
}
//=====
void TPCDigitsContainer::Reset()
{
  // clears counters and resets buffer's entries
  for(int i=0; i!=72; ++i) {
    for(unsigned int j=0; j!=fDigits[i].size(); ++j) (fDigits[i].at(j))->Reset();
    fNDigits[i] = 0;
  }
}
//=====
void TPCDigitsContainer::Add(Module_t mod, Pad_t pad, BinTime_t tim, Adc_t w, int tid)
{
  // adds digit into container
  TPCDigit *mine;
  // attempt to find the pad in buffer
  std::map<Pad_t,TPCDigit*>::const_iterator iter = fDigits[mod].find(pad);
  if(iter==fDigits[mod].end()) {
    // not found. new digit must be added
    mine = new TPCDigit();
    (fDigits[mod])[pad] = mine; // adding new key
    fNDigits[mod]++;
  } else {
    // found. reuse digit and increase
    // inc counter only if it was empty
    mine = (*iter).second;
    if(mine->GetEntries()==0) fNDigits[mod]++;
  }
  // filling new weight
  mine->AddContent(tim,w,tid);
}
//=====
TPCDigit* TPCDigitsContainer::GetDigit(Module_t mod, Pad_t pad)
{
  // return either digit or NULL
  std::map<Pad_t,TPCDigit*>::const_iterator iter = fDigits[mod].find(pad);
  if(iter==fDigits[mod].end()) return NULL;
  TPCDigit *me = (*iter).second;
  if(me->GetEntries()==0) return NULL;
  return me;
}
