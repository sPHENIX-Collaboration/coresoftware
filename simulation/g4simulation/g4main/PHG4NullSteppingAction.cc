#include "PHG4NullSteppingAction.h"
#include "PHG4EmcDetector.h"

#include <PHG4HitContainer.h>
#include <PHG4Hit.h>
#include <PHG4Hitv1.h>
#include <G4Step.hh>

#include <getClass.h>

#include <iostream>

using namespace std;

//____________________________________________________________________________..
bool 
PHG4NullSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  //G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  //std::cout << "PHG4NullSteppingAction::UserSteppingAction: Stepping in volume " << volume->GetName() << std::endl;
  return false;
}
