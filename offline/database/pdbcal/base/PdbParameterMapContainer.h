#ifndef PDBCAL_BASE_PDBPARAMETERMAPCONTAINER_H
#define PDBCAL_BASE_PDBPARAMETERMAPCONTAINER_H

#include "PdbCalChan.h"

#include <map>
#include <string>
#include <utility>

class PdbParameterMap;

class PdbParameterMapContainer : public PdbCalChan
{
 public:
  typedef std::map<int, PdbParameterMap *> parMap;
  typedef parMap::const_iterator parIter;
  typedef std::pair<parIter, parIter> parConstRange;

  PdbParameterMapContainer() {}
  ~PdbParameterMapContainer() override;

  void print() const override;

  void Reset() override;

  void AddPdbParameterMap(const int layer, PdbParameterMap *params);
  const PdbParameterMap *GetParameters(const int layer) const;
  PdbParameterMap *GetParametersToModify(const int layer);
  parConstRange get_ParameterMaps() const { return make_pair(parametermap.begin(), parametermap.end()); }

  //! write PdbParameterMapContainer to an external file with root or xml extension.
  int WriteToFile(const std::string &detector_name, const std::string &extension, const std::string &dir = ".");

 protected:
  std::map<int, PdbParameterMap *> parametermap;

  ClassDefOverride(PdbParameterMapContainer, 1)
};

#endif
