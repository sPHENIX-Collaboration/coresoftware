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

  virtual void identify(std::ostream& = std::cout) const override;
  virtual std::size_t size() const;

  int LoadFromFile(std::string const& = "InttSurveyMap.root");
  int LoadFromCDB(std::string const& = "InttSurveyMap");

  virtual val_t const* GetAbsoluteTransform(key_t) const;
  virtual val_t const* GetRelativeTransform(key_t) const;

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

  virtual int v_LookupAbsoluteTransform(key_t const&, map_t::const_iterator&) const;
  virtual int v_LookupRelativeTransform(key_t const&, map_t::const_iterator&) const;

 private:
  ClassDefOverride(InttSurveyMap, 1)
};

#endif  // INTT_SURVEY_MAP_H
