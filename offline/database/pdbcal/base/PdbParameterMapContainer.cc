#include "PdbParameterMapContainer.h"
#include "PdbParameterMap.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <iostream>

using namespace std;

PdbParameterMapContainer::~PdbParameterMapContainer()
{
  while(parametermap.begin() != parametermap.end())
    {
      delete parametermap.begin()->second;
      parametermap.erase(parametermap.begin());
    }
  return;
}

void
PdbParameterMapContainer::print() const
{
  for (map<int, PdbParameterMap *>::const_iterator iter = parametermap.begin(); 
       iter != parametermap.end(); ++iter)
    {
      cout << "layer " << iter->first << endl;
      iter->second->print();
    }
  return;
}

void
PdbParameterMapContainer::AddPdbParameterMap(const int layer, PdbParameterMap *params)
{
  if (parametermap.find(layer) != parametermap.end())
    {
      cout << PHWHERE << " layer " << layer << " already exists" << endl; 
      gSystem->Exit(1);
    }
  parametermap[layer] = params;
}

const PdbParameterMap *
PdbParameterMapContainer::GetParameters(const int layer) const
{
  map<int, PdbParameterMap *>::const_iterator iter = parametermap.find(layer);
if (iter == parametermap.end())
  {
    return NULL;
  }
 return iter->second;
}

PdbParameterMap *
PdbParameterMapContainer::GetParametersToModify(const int layer)
{
  map<int, PdbParameterMap *>::iterator iter = parametermap.find(layer);
  if (iter == parametermap.end())
    {
      return NULL;
    }
  return iter->second;
}
