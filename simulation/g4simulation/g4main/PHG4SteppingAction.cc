/*!
 * \file PHG4SteppingAction.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2015/01/07 23:50:05 $
 */

#include "PHG4SteppingAction.h"

#include <Geant4/G4Step.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4EmSaturation.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>
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
