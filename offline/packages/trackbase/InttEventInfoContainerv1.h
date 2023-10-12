#ifndef INTT_EVENT_INFO_CONTAINERv1_H
#define INTT_EVENT_INFO_CONTAINERv1_H

#include "InttEventInfoContainer.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>
#include <string>
#include <map>

class InttEventInfoContainerv1 : public InttEventInfoContainer 
{
public:
  InttEventInfoContainerv1();
  ~InttEventInfoContainerv1() override;

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override;

  void AddInfo(KEY_t const&, VAL_t const&) override;
  VAL_t& GetInfo(KEY_t const&) override;

protected:
  Map info_map;

 private:
  ClassDefOverride(InttEventInfoContainer, 1)
};

#endif//INTT_EVENT_INFO_CONTAINERv1_H
