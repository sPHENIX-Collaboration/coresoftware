/*!
 * \file PHG4SteppingAction.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2015/01/07 23:50:05 $
 */

#include "PHG4SteppingAction.h"
#include "PHG4Hit.h"

#include <Geant4/G4Step.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4EmSaturation.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4StepPoint.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4NavigationHistory.hh>



#include <iostream>
#include <cassert>

using namespace std;

double
PHG4SteppingAction::GetScintLightYield(const G4Step* step)
{
  //      pirated from G4Scintillation::PostStepDoIt()

  double light_yield = 0;

  const G4Track* aTrack = step->GetTrack();
  const G4Material* aMaterial = aTrack->GetMaterial();
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
      aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable)
    {
      string mname(aMaterial->GetName());

      std::set<std::string>::const_iterator it =
          _ScintLightYieldMissingMaterial.find(mname);

      if (it == _ScintLightYieldMissingMaterial.end())
        {
          _ScintLightYieldMissingMaterial.insert(mname);

          cout << "PHG4SteppingAction::GetScintLightYield - WARNING - "
              << "can not find Material Properties Table for material " << mname
              << ", will assume it do NOT scintillate. "
              << "Please ignore this warning if you do not expect scintillation light from "
              << mname << endl;
        }

      return 0;
    }

  if (aMaterialPropertiesTable->ConstPropertyExists("SCINTILLATIONYIELD"))
    {
      light_yield = aMaterialPropertiesTable->GetConstProperty(
          "SCINTILLATIONYIELD") * GetVisibleEnergyDeposition(step) * GeV;

      return light_yield;
    } //  if (aMaterialPropertiesTable->ConstPropertyExists("SCINTILLATIONYIELD"))
  else
    {

      string mname(aMaterial->GetName());

      std::set<std::string>::const_iterator it =
          _ScintLightYieldMissingMaterial.find(mname);

      if (it == _ScintLightYieldMissingMaterial.end())
        {
          _ScintLightYieldMissingMaterial.insert(mname);

          cout << "PHG4SteppingAction::GetScintLightYield - WARNING - "
              << "can not find scintillation light yield for material " << mname
              << ", will assume it do NOT scintillate. "
              << "Please ignore this warning if you do not expect scintillation light from "
              << mname << endl;
        }

      return 0;
    }

  return light_yield;
}

double
PHG4SteppingAction::GetVisibleEnergyDeposition(const G4Step* step)
{

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  if (emSaturation)
    {
      if (verbosity)
	{
	  emSaturation->SetVerbose(verbosity);
	}
      double visen = emSaturation->VisibleEnergyDepositionAtAStep(step) / GeV;
      return visen;
    }
  else
    {
      cout
          << "PHG4SteppingAction::GetScintLightYield - ERROR - can NOT initialize G4EmSaturation!"
          << endl;

      return 0;
    }
}

void
PHG4SteppingAction::StoreLocalCoordinate(PHG4Hit * hit, const G4Step* aStep,
    const bool do_prepoint, const bool do_postpoint)
{
  assert(hit);
  assert(aStep);

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();

  if (do_prepoint)
    {

      G4ThreeVector worldPosition = preStepPoint->GetPosition();
      G4ThreeVector localPosition =
          theTouchable->GetHistory()->GetTopTransform().TransformPoint(
              worldPosition);

      hit->set_local_x(0, localPosition.x() / cm);
      hit->set_local_y(0, localPosition.y() / cm);
      hit->set_local_z(0, localPosition.z() / cm);
    }
  if (do_postpoint)
    {
      G4StepPoint * postPoint = aStep->GetPostStepPoint();

      G4ThreeVector worldPosition = postPoint->GetPosition();
      G4ThreeVector localPosition =
          theTouchable->GetHistory()->GetTopTransform().TransformPoint(
              worldPosition);

      hit->set_local_x(1, localPosition.x() / cm);
      hit->set_local_y(1, localPosition.y() / cm);
      hit->set_local_z(1, localPosition.z() / cm);
    }

}

bool
PHG4SteppingAction::IntOptExist(const std::string &name)
{
  if (opt_int.find(name) != opt_int.end())
    {
      return true;
    }
  return false;
}

int
PHG4SteppingAction::GetIntOpt(const std::string &name)
{
  if (IntOptExist(name))
    {
      return opt_int.find(name)->second;
    }
  cout << "option " << name << " does not exist" << endl;
  exit(1);
}
