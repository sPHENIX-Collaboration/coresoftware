#include "PHG4StepStatusDecode.h"

#include <Geant4/G4StepStatus.hh>

#include <TSystem.h>

#include <iostream>
#include <map>

static std::map<int, std::string> stepstatus;

std::string PHG4StepStatusDecode::GetStepStatus(const int istatus)
{
  if (stepstatus.empty())
  {
    //    std::cout << "filling stepstatus map" << std::endl;
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
    std::cout << "could not find status " << istatus << " in stepstatus map" << std::endl;
    gSystem->Exit(1);
  }
  return stepstatus[istatus];
}
