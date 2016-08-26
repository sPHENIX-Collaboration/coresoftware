#include "PHG4ParametersContainer.h"
#include "PHG4Parameters.h"

#include <phool/phool.h>

#include <TSystem.h>

#include <iostream>

using namespace std;

PHG4ParametersContainer::~PHG4ParametersContainer()
{
  while(parametermap.begin() != parametermap.end())
    {
      delete parametermap.begin()->second;
      parametermap.erase(parametermap.begin());
    }

}

void
PHG4ParametersContainer::AddPHG4Parameters(const int layer, PHG4Parameters *params)
{
  if (parametermap.find(layer) != parametermap.end())
    {
      cout << PHWHERE << " layer " << layer << " already exists for " 
	   << (parametermap.find(layer))->second->Name() << endl;
      gSystem->Exit(1);
    }
  parametermap[layer] = params;
}

const PHG4Parameters *
PHG4ParametersContainer::GetParameters(const int layer) const
{
  map<int, PHG4Parameters *>::const_iterator iter = parametermap.find(layer);
if (iter == parametermap.end())
  {
    cout << "could not find parameters for layer " << layer
	 << endl;
    return NULL;
  }
 return iter->second;
}

PHG4Parameters *
PHG4ParametersContainer::GetParametersToModify(const int layer)
{
  map<int, PHG4Parameters *>::iterator iter = parametermap.find(layer);
  if (iter == parametermap.end())
    {
      cout << "could not find parameters for layer " << layer
	   << endl;
      return NULL;
    }
  return iter->second;
}
