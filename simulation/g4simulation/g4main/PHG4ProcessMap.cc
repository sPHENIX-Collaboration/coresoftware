
#include "PHG4ProcessMap.h"
#include "PHG4MCProcessDefs.h"

#include <Geant4/G4VProcess.hh>

#include <iomanip>
#include <iostream>
#include <string>

//
// private methods
//

//_____________________________________________________________________________
bool PHG4ProcessMap::IsDefined(int subType)
{
  return !(fMap.find(subType) == fMap.end());
}

//
// public methods
//

//_____________________________________________________________________________
bool PHG4ProcessMap::Add(int subType, PHG4MCProcess mcProcess)
{
  if (!IsDefined(subType))
  {
    // insert into map
    // only in case it is not yet here
    fMap[subType] = mcProcess;
    return true;
  }
  return false;
}

//_____________________________________________________________________________
void PHG4ProcessMap::PrintAll() const
{
  if (fMap.empty())
  {
    return;
  }

  std::cout << "Dump of PHG4ProcessMap - " << fMap.size()
            << " entries:" << std::endl;
  int counter = 0;
  for (auto [subType, codes] : fMap)
  {
    // TO DO: get process sub-type name
    std::cout << "Map element " << std::setw(3) << counter++ << "   "
              << subType << "   " << PHG4MCProcessName[codes] << std::endl;
  }
}

//_____________________________________________________________________________
void PHG4ProcessMap::Clear()
{
  if (fMap.size())
  {
    fMap.clear();
  }
}

//_____________________________________________________________________________
PHG4MCProcess
PHG4ProcessMap::GetMCProcess(const G4VProcess* process) const
{
  if (!process)
  {
    return kPNoProcess;
  }

  auto i = fMap.find(process->GetProcessSubType());
  if (i == fMap.end())
  {
    std::string text = "Unknown process code for ";
    text += process->GetProcessName();
    std::cerr << "PHG4ProcessMap::GetCodes " << text.c_str() << std::endl;
    return kPNoProcess;
  }
  else
  {
    return (*i).second;
  }
}

//_____________________________________________________________________________
std::string_view PHG4ProcessMap::GetMCProcessName(const G4VProcess* process) const
{
  if (!process)
  {
    return PHG4MCProcessName[kPNoProcess];
  }

  return PHG4MCProcessName[GetMCProcess(process)];
}
