#ifndef INTT_FEE_MAPv1_H
#define INTT_FEE_MAPv1_H

#include "InttMap.h"
#include "InttFeeMap.h"

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <string>

class CDBTTree;

class InttFeeMapv1 : public InttFeeMap
{
 public:
  InttFeeMapv1() = default;
  ~InttFeeMapv1() override;

  void identify(std::ostream& = std::cout) const override;

  int Convert(InttMap::Online_s&, InttMap::Offline_s const&) const override;
  int Convert(InttMap::Offline_s&, InttMap::Online_s const&) const override;

  int Convert(InttMap::RawData_s&, InttMap::Online_s const&) const override;
  int Convert(InttMap::Online_s&, InttMap::RawData_s const&) const override;

  int Convert(InttMap::RawData_s&, InttMap::Offline_s const&) const override;
  int Convert(InttMap::Offline_s&, InttMap::RawData_s const&) const override;

 protected:
  int v_LoadFromCDBTTree(CDBTTree&) override;

 private:
  std::map<InttMap::Online_s, InttMap::RawData_s, InttMap::OnlineWildcardComparator>* m_onl_to_raw{nullptr};
  std::map<InttMap::RawData_s, InttMap::Online_s, InttMap::RawDataWildcardComparator>* m_raw_to_onl{nullptr};

  ClassDefOverride(InttFeeMapv1, 1)
};

#endif  // INTT_FEE_MAPv1_H
