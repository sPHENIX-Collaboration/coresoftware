#ifndef PdbParameterMapContainer__h
#define PdbParameterMapContainer__h

#include "PdbCalChan.h"

#include <map>

class PdbParameterMap;

class PdbParameterMapContainer: public PdbCalChan
{
 public:

  typedef std::map<int, PdbParameterMap *> parMap;
  typedef parMap::const_iterator parIter;
  typedef std::pair<parIter, parIter> parConstRange;

  PdbParameterMapContainer() {}
  virtual ~PdbParameterMapContainer();

  void print() const;

  void AddPdbParameterMap(const int layer, PdbParameterMap *params);
  const PdbParameterMap *GetParameters(const int layer) const;
  PdbParameterMap *GetParametersToModify(const int layer);
  parConstRange get_ParameterMaps() const {return make_pair(parametermap.begin(), parametermap.end());}

 protected:
  std::map<int, PdbParameterMap *> parametermap;

  ClassDef(PdbParameterMapContainer,1)
};

#endif
