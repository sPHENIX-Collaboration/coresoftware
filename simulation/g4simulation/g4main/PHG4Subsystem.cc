#include "PHG4Subsystem.h"

#include <TSystem.h>

#include <iostream>

using namespace std;

void PHG4Subsystem::SetMotherSubsystem(PHG4Subsystem *subsys)
{
  if (subsys->CanBeMotherSubsystem())
  {
     m_MyMotherSubsystem = subsys;
     return;
  }
  cout << "PHG4Subsystem::SetMotherSubsystem: "
       << subsys->Name() << " is not implemented as a mother subsystem"
       << endl;
  gSystem->Exit(1);
}
