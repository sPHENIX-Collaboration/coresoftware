/*!
 * \file PHG4SteppingAction.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.1 $
 * \date $Date: 2015/01/07 23:50:05 $
 */

#include "PHG4SteppingAction.h"
#include "PHG4Hit.h"

#include <Geant4/G4EmSaturation.hh>
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4Track.hh>

#include <algorithm>
#include <cassert>
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
        m_ScintLightYieldMissingMaterialSet.find(mname);

    if (it == m_ScintLightYieldMissingMaterialSet.end())
    {
      m_ScintLightYieldMissingMaterialSet.insert(mname);

      cout << "PHG4SteppingAction::GetScintLightYield - WARNING - "
           << "can not find Material Properties Table for material " << mname
           << ", will assume it do NOT scintillate. "
           << "Please ignore this warning if you do not expect scintillation light from "
           << mname << endl;
    }

    return 0.;
  }

  if (aMaterialPropertiesTable->ConstPropertyExists("SCINTILLATIONYIELD"))
  {
    light_yield = aMaterialPropertiesTable->GetConstProperty(
                      "SCINTILLATIONYIELD") *
                  GetVisibleEnergyDeposition(step) * GeV;

    return light_yield;
  }  //  if (aMaterialPropertiesTable->ConstPropertyExists("SCINTILLATIONYIELD"))
  else
  {
    string mname(aMaterial->GetName());

    std::set<std::string>::const_iterator it =
        m_ScintLightYieldMissingMaterialSet.find(mname);

    if (it == m_ScintLightYieldMissingMaterialSet.end())
    {
      m_ScintLightYieldMissingMaterialSet.insert(mname);

      cout << "PHG4SteppingAction::GetScintLightYield - WARNING - "
           << "can not find scintillation light yield for material " << mname
           << ", will assume it do NOT scintillate. "
           << "Please ignore this warning if you do not expect scintillation light from "
           << mname << endl;
    }

    return 0.;
  }

  return light_yield;
}

double
PHG4SteppingAction::GetVisibleEnergyDeposition(const G4Step* step)
{
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  if (emSaturation)
  {
    if (m_Verbosity)
    {
      emSaturation->SetVerbose(m_Verbosity);
    }
    double visen = emSaturation->VisibleEnergyDepositionAtAStep(step) / GeV;
    return visen;
  }
  else
  {
    cout
        << "PHG4SteppingAction::GetScintLightYield - ERROR - can NOT initialize G4EmSaturation!"
        << endl;

    return 0.;
  }
}

void PHG4SteppingAction::StoreLocalCoordinate(PHG4Hit* hit, const G4Step* aStep,
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
    G4StepPoint* postPoint = aStep->GetPostStepPoint();

    G4ThreeVector worldPosition = postPoint->GetPosition();
    G4ThreeVector localPosition =
        theTouchable->GetHistory()->GetTopTransform().TransformPoint(
            worldPosition);

    hit->set_local_x(1, localPosition.x() / cm);
    hit->set_local_y(1, localPosition.y() / cm);
    hit->set_local_z(1, localPosition.z() / cm);
  }
}

void PHG4SteppingAction::SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr)
{
  m_LightBalanceInnerRadius = inner_radius;
  m_LightBalanceInnerCorr = inner_corr;
  m_LightBalanceOuterRadius = outer_radius;
  m_LightBalanceOuterCorr = outer_corr;
  return;
}

double PHG4SteppingAction::GetLightCorrection(const double xpos, const double ypos) const
{
  double r = sqrt(xpos * xpos + ypos * ypos);
  double correction = GetLightCorrection(r);
  return correction;
}

double PHG4SteppingAction::GetLightCorrection(const double r) const
{
  double correction = 1.;
  if (ValidCorrection())
  {
    double m = (m_LightBalanceOuterCorr - m_LightBalanceInnerCorr) / (m_LightBalanceOuterRadius - m_LightBalanceInnerRadius);
    double b = m_LightBalanceInnerCorr - m * m_LightBalanceInnerRadius;
    correction = m * r + b;
    correction = std::min(correction, 1.);
    correction = std::max(correction, 0.);
  }
  return correction;
}

bool PHG4SteppingAction::ValidCorrection() const
{
  if (isfinite(m_LightBalanceOuterRadius) &&
      isfinite(m_LightBalanceInnerRadius) &&
      isfinite(m_LightBalanceOuterCorr) &&
      isfinite(m_LightBalanceInnerCorr))
  {
    return true;
  }
  return false;
}
