#ifndef PdbParameterMapContainer__h
#define PdbParameterMapContainer__h

#include "PdbCalChan.h"

#include <map>

class PdbParameterMap;

class PdbParameterMapContainer: public PdbCalChan
{
 public:
  PdbParameterMapContainer() {}
  virtual ~PdbParameterMapContainer();

  void print() const;

  void AddPdbParameterMap(const int layer, PdbParameterMap *params);
  const PdbParameterMap *GetParameters(const int layer) const;
  PdbParameterMap *GetParametersToModify(const int layer);

 protected:
  std::map<int, PdbParameterMap *> parametermap;

  ClassDef(PdbParameterMapContainer,1)
};

#endif
