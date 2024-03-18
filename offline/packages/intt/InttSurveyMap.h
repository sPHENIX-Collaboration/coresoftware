#ifndef INTT_SURVEY_MAP_H
#define INTT_SURVEY_MAP_H

#include "InttMap.h"

#include <phool/PHObject.h>

#include <Eigen/Geometry>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <string>

class CDBTTree;

class InttSurveyMap : public PHObject
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

  virtual void identify(std::ostream& = std::cout) const override;
  virtual std::size_t size() const;

  virtual val_t const* GetAbsoluteTransform(key_t const&) const;
  virtual val_t const* GetRelativeTransform(key_t const&) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  ClassDefOverride(InttSurveyMap, 1)
};

#endif  // INTT_SURVEY_MAP_H
