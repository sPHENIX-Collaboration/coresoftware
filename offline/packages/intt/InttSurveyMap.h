#ifndef INTT_SURVEY_MAP_H
#define INTT_SURVEY_MAP_H

#include "InttLoadable.h"
#include "InttMap.h"

#include <Eigen/Geometry>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <string>

class CDBTTree;

class InttSurveyMap : public InttLoadable
{
 public:
  typedef std::map<InttMap::Offline_s, Eigen::Affine3d, InttMap::OfflineComparator> map_t;
  typedef InttMap::Offline_s key_t;
  typedef Eigen::Affine3d val_t;

  InttSurveyMap() = default;
  ~InttSurveyMap() override = default;

  int LoadFromFile(std::string const& = "InttSurveyMap.root");
  int LoadFromCDB(std::string const& = "InttSurveyMap");

  int GetStripTransform(key_t const&, val_t&) const;
  int GetSensorTransform(key_t const&, val_t&) const;
  int GetLadderTransform(key_t const&, val_t&) const;

  virtual val_t const* GetAbsoluteTransform(key_t const&) const;

 protected:
  int LoadFromCdbTTree(CDBTTree&) override;

 private:
  map_t m_absolute_transforms;
};

#endif  // INTT_SURVEY_MAP_H
