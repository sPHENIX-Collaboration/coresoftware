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

  val_t const* GetAbsoluteTransform(key_t const&) const override;

 protected:
  int v_LoadFromCDBTTree(CDBTTree&) override;

 private:
  map_t* m_absolute_transforms = nullptr;

  ClassDefOverride(InttSurveyMapv1, 1)
};

#endif  // INTT_SURVEY_MAPv1_H
