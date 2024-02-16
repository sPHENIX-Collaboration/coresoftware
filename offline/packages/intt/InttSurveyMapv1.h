#ifndef INTT_SURVEY_MAPv1_H
#define INTT_SURVEY_MAPv1_H

#include "InttSurveyMap.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>  // for map<>::const_iterator

class InttSurveyMapv1 : public InttSurveyMap
{
 public:
  using InttSurveyMap::key_t;
  using InttSurveyMap::map_t;
  using InttSurveyMap::val_t;

  InttSurveyMapv1() = default;
  ~InttSurveyMapv1() override;

  void identify(std::ostream& = std::cout) const override;
  std::size_t size() const override;

 protected:
  int v_LoadFromCDBTTree(CDBTTree&) override;

  int v_LookupAbsoluteTransform(key_t const&, map_t::const_iterator&) const override;
  int v_LookupRelativeTransform(key_t const&, map_t::const_iterator&) const override;

 private:
  map_t* m_absolute_transforms = nullptr;
  map_t* m_relative_transforms = nullptr;

  ClassDefOverride(InttSurveyMapv1, 1)
};

#endif  // INTT_SURVEY_MAPv1_H
