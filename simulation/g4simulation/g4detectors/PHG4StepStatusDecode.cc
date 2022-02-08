#include "PHG4StepStatusDecode.h"

#include <Geant4/G4StepStatus.hh>

#include <TSystem.h>

#include <iostream>
#include <map>

using namespace std;
static map<int, string> stepstatus;

string PHG4StepStatusDecode::GetStepStatus(const int istatus)
{
  if (stepstatus.empty())
  {
    cout << "filling stepstatus map" << endl;
    stepstatus[fWorldBoundary] = "fWorldBoundary";
    stepstatus[fGeomBoundary] = "fGeomBoundary";
    stepstatus[fAtRestDoItProc] = "fAtRestDoItProc";
    stepstatus[fAlongStepDoItProc] = "fAlongStepDoItProc";
    stepstatus[fPostStepDoItProc] = "fPostStepDoItProc";
    stepstatus[fUserDefinedLimit] = "fUserDefinedLimit";
    stepstatus[fExclusivelyForcedProc] = "fExclusivelyForcedProc";
    stepstatus[fUndefined] = "fUndefined";
  }
  if (stepstatus.find(istatus) == stepstatus.end())
  {
    cout << "could not find status " << istatus << " in stepstatus map" << endl;
    gSystem->Exit(1);
  }
  return stepstatus[istatus];
}
