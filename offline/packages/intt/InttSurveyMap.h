#ifndef INTT_SURVEY_MAP_H
#define INTT_SURVEY_MAP_H

#include "InttMap.h"
#include "InttLoadable.h"

#include <Eigen/Geometry>

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
  virtual ~InttSurveyMap() override;

  int GetStripTransform(key_t, val_t&) const;
  int GetSensorTransform(key_t, val_t&) const;
  int GetLadderTransform(key_t, val_t&) const;

  val_t const* GetTransform(key_t const&) const;

  std::string DefaultName() const override { return "InttSurveyMap"; }

 protected:
  int LoadFromCDBTTree(CDBTTree&) override;

 private:
  map_t* m_transforms = nullptr;
};

#endif  // INTT_SURVEY_MAP_H
