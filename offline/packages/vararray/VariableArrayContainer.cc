#include "VariableArrayContainer.h"

#include "VariableArray.h"

#include <TObjArray.h>

using namespace std;

VariableArrayContainer::VariableArrayContainer()
{
  arraycontainer = new TObjArray();
  return;
}

VariableArrayContainer::~VariableArrayContainer()
{
  TObjArrayIter *iter = new TObjArrayIter(arraycontainer);
  while (VariableArray *vararray = dynamic_cast<VariableArray *>(iter->Next()))
  {
    delete vararray;
  }
  delete iter;
  delete arraycontainer;
  return;
}

void VariableArrayContainer::identify(ostream &os) const
{
  os << "contain ObjCont" << endl;
  return;
}

void VariableArrayContainer::AddVarArray(VariableArray *var)
{
  arraycontainer->Add(var);
  cout << "Adding " << var->GetName() << endl;
  return;
}

void VariableArrayContainer::Reset()
{
  TObjArrayIter *iter = new TObjArrayIter(arraycontainer);
  while (VariableArray *vararray = dynamic_cast<VariableArray *>(iter->Next()))
  {
    vararray->Reset();
  }
  delete iter;
  return;
}
