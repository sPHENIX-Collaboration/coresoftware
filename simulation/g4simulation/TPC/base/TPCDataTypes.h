// TPC definitions
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
  typedef UChar_t Module_t;
  typedef UShort_t Pad_t;
  typedef UShort_t Time_t;
  typedef UShort_t Adc_t;
  typedef std::pair<Pad_t,Float_t> PadQuota_t;
  typedef std::pair<Time_t,Adc_t> TimeQuota_t;
  typedef std::vector<Pad_t> PadRange_t;
  typedef std::vector<Module_t> ModuleRange_t;
  typedef std::vector<PadQuota_t> PadQuotaRange_t;
  typedef std::vector<TimeQuota_t> TimeQuotaRange_t;
  typedef std::pair<Float_t,Float_t> PairOfFloats_t;
  typedef std::map<Pad_t,PairOfFloats_t> MapPadXY_t;
};

#endif
