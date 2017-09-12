// holds the data types used for the TPC
// Author: Carlos Perez
#ifndef __TPCDATATYPES_H__
#define __TPCDATATYPES_H__

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include <vector>
#include <map>
#include <utility>

namespace TPCDataTypes {
  // basal datatypes
  typedef unsigned char  Module_t;
  typedef unsigned short Pad_t;
  typedef unsigned short BinTime_t;
  typedef unsigned short Adc_t;

  // derivated datatypes
  typedef std::pair<Pad_t,Float_t> PadQuota_t;
  typedef std::pair<BinTime_t,Adc_t> TimeQuota_t;
  typedef std::vector<Pad_t> PadRange_t;
  typedef std::vector<Module_t> ModuleRange_t;
  typedef std::vector<PadQuota_t> PadQuotaRange_t;
  typedef std::vector<TimeQuota_t> TimeQuotaRange_t;
  typedef std::pair<Float_t,Float_t> PairOfFloats_t;
  typedef std::map<Pad_t,PairOfFloats_t> MapPadXY_t;
  typedef std::pair<Adc_t,std::map<int,int>> AdcTid_t;

  typedef std::map<BinTime_t,AdcTid_t> TrainOfDigits_t;
  typedef std::map<BinTime_t,AdcTid_t>::iterator TrainOfDigits_Iter_t;
  typedef std::map<BinTime_t,AdcTid_t>::const_iterator TrainOfDigits_ConstIter_t;
};

#endif
