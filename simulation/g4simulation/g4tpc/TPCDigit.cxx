// TPC DIGIT class
// Stores vector of digits associated to one readout pad
// Author: Carlos Perez
#include <map>
#include "TPCDigit.h"
#include "TPCDataTypes.h"

using namespace TPCDataTypes;

//=====
TPCDigit::TPCDigit() {}
//=====
TPCDigit::~TPCDigit() {
  TrainOfDigits_t *tod = fCell->get_train_of_digits();
  tod->clear();
}
//=====
Adc_t TPCDigit::GetContent(const BinTime_t bin) {
  TrainOfDigits_t *tod = fCell->get_train_of_digits();
  TrainOfDigits_ConstIter_t iter = tod->find(bin);
  if(iter==tod->end()) return 0;
  return ((*iter).second).first;
}
//=====
void TPCDigit::AddContent(const BinTime_t bin, Adc_t val, int tid) {
  TrainOfDigits_t *tod = fCell->get_train_of_digits();
  TrainOfDigits_Iter_t iter = tod->find(bin);
  if(iter==tod->end()) {
    // create new
    AdcTid_t entry;
    entry.first = val;
    entry.second[tid] = 1;
    (*tod)[bin] = entry;
  } else {
    // add to existing
    ((*iter).second).first += val;
    ((*iter).second).second[tid] += 1;
  }
}
//=====
int TPCDigit::GetEntries() {
  TrainOfDigits_t *tod = fCell->get_train_of_digits();
  return static_cast<int>(tod->size());
}
//=====
void TPCDigit::Reset() {
  TrainOfDigits_t *tod = fCell->get_train_of_digits();
  tod->clear();
}
