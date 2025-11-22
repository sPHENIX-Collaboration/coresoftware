#include "PHG4Subsystem.h"

#include <TSystem.h>

#include <iostream>

void PHG4Subsystem::SetMotherSubsystem(PHG4Subsystem *subsys)
{
  if (subsys->CanBeMotherSubsystem())
  {
    m_MyMotherSubsystem = subsys;
    return;
  }
  std::cout << "PHG4Subsystem::SetMotherSubsystem: "
            << subsys->Name() << " is not implemented as a mother subsystem"
            << std::endl;
  gSystem->Exit(1);
}
