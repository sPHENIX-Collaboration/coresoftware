#include "PHG4ZDCSteppingAction.h"

#include "PHG4ZDCDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <Geant4/G4DynamicParticle.hh>  // for G4DynamicParticle
#include <Geant4/G4IonisParamMat.hh>    // for G4IonisParamMat
#include <Geant4/G4Material.hh>         // for G4Material
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>  // for basic_string, operator+

class PHCompositeNode;

//____________________________________________________________________________..
PHG4ZDCSteppingAction::PHG4ZDCSteppingAction(PHG4ZDCDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_IsActiveFlag(m_Params->get_int_param("active"))
  , absorbertruth(m_Params->get_int_param("absorberactive"))
  , m_IsBlackHole(m_Params->get_int_param("blackhole"))
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int seed = PHRandomSeed();
  gsl_rng_set(RandomGenerator, seed);
}

PHG4ZDCSteppingAction::~PHG4ZDCSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
  gsl_rng_free(RandomGenerator);
}

//____________________________________________________________________________..
bool PHG4ZDCSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInZDC(volume)
  // returns
  //  0 is outside of ZDC
  //  1 is inside scintillator
  // -1 is inside absorber or support structure (dead material)

  int whichactive = m_Detector->IsInZDC(volume);

  if (!whichactive)
  {
    return false;
  }

  int layer_id = m_Detector->get_Layer();
  int idx_k = -1;
  int idx_j = -1;

  if (whichactive > 0)  // in scintillator or fiber
  {
    /* Find indices of scintillator / tower containing this step */
    //FindIndex(touch, idx_j, idx_k);
    if (whichactive == 2) FindIndexZDC(touch, idx_j, idx_k);
    if (whichactive == 1) FindIndexSMD(touch, idx_j, idx_k);
  }
  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (m_IsBlackHole)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  /* Make sure we are in a volume */
  if (m_IsActiveFlag)
  {
    /* Check if particle is 'geantino' */
    bool geantino = false;
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
        aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
    {
      geantino = true;
    }

    /* Get Geant4 pre- and post-step points */
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();

    //if prepoint is in fiber
    if (whichactive > 1)
    {
      double charge = aTrack->GetParticleDefinition()->GetPDGCharge();
      //if charged particle
      if (charge != 0)
      {
        //check if prepoint in active volume & postpoint out of active volume
        G4VPhysicalVolume* postvolume = postPoint->GetTouchableHandle()->GetVolume();
        int postactive = m_Detector->IsInZDC(postvolume);
        //postpoint outside fiber
        if (!postactive)
        {
          //get particle information here
          int pid = aTrack->GetParticleDefinition()->GetPDGEncoding();
          //calculate incidence angle
          const G4DynamicParticle* dypar = aTrack->GetDynamicParticle();
          const G4ThreeVector& pdirect = dypar->GetMomentumDirection();
          // this triggers cppcheck, the code is good and the warning is suppressed
          // cppcheck-suppress [duplicateAssignExpression, unmatchedSuppression]
          double dy = sqrt(2) / 2.;
          double dz = sqrt(2) / 2.;
          if (idx_j == 1) dz = -dz;
          double CosTheta = pdirect.y() * dy + pdirect.z() * dz;
          double angle = acos(CosTheta) * 180.0 / M_PI;
          if (pid == 11 || pid == -11)
          {
            //find energy
            G4double E = dypar->GetTotalEnergy();
            //electron response here
            double avg_ph = ZDCEResponse(E, angle);
            avg_ph *= 0.16848;
            //use Poisson Distribution here
            int n_ph = gsl_ran_poisson(RandomGenerator, avg_ph);
            light_yield += n_ph;
          }
          else
          {
            G4double E = dypar->GetTotalEnergy();
            G4double P = dypar->GetTotalMomentum();
            double beta = P / E;
            double avg_ph = ZDCResponse(beta, angle);
            avg_ph *= 0.16848;
            int n_ph = gsl_ran_poisson(RandomGenerator, avg_ph);
            light_yield += n_ph;
          }
        }
      }
    }
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }

      /* Set hit location (space point) */
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      /* Set hit time */
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

      //set the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      /* set intial energy deposit */
      m_Hit->set_edep(0);

      /* Now add the hit to the hit collection */
      // here we do things which are different between scintillator and absorber hits
      if (whichactive > 0)
      {
        m_CurrentHitContainer = m_HitContainer;
        m_Hit->set_eion(0);
        m_Hit->set_light_yield(0);  // for scintillator only, initialize light yields
        /* Set hit location (tower index) */
        m_Hit->set_index_k(idx_k);
        m_Hit->set_index_j(idx_j);
      }
      else
      {
        if (whichactive == -1)
        {
          m_CurrentHitContainer = m_AbsorberHitContainer;
        }
        else
        {
          m_CurrentHitContainer = m_SupportHitContainer;
        }
      }
      // here we set what is common for scintillator and absorber hits
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_CurrentShower = pp->GetShower();
        }
      }
      break;
    default:
      break;
    }

    if (whichactive == 1)
    {
      //light_yield = eion;
      light_yield = GetVisibleEnergyDeposition(aStep);  // for scintillator only, calculate light yields
      static bool once = true;

      if (once && edep > 0)
      {
        once = false;

        if (Verbosity() > 0)
        {
          std::cout << "PHG4ZDCSteppingAction::UserSteppingAction::"
                    //
                    << m_Detector->GetName() << " - "
                    << " use scintillating light model at each Geant4 steps. "
                    << "First step: "
                    << "Material = "
                    << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetName()
                    << ", "
                    << "Birk Constant = "
                    << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetIonisation()->GetBirksConstant()
                    << ","
                    << "edep = " << edep << ", "
                    << "eion = " << eion
                    << ", "
                    << "light_yield = " << light_yield << std::endl;
        }
      }
    }

    /* sum up the energy to get total deposited */

    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (whichactive > 0)
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    }

    // if any of these conditions is true this is the last step in
    // this volume and we need to save the hit
    // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
    // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
    // (not sure if this will ever be the case)
    // aTrack->GetTrackStatus() == fStopAndKill: track ends
    if (postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // Update exit values
      m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
      m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
      m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

      m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

      // special case for geantinos
      if (geantino)
      {
        m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
        if (whichactive > 0)
        {
          m_Hit->set_eion(-1);
          m_Hit->set_light_yield(-1);
        }
      }
      // save only hits with energy deposit (or -1 for geantino)
      if (m_Hit->get_edep())
      {
        m_CurrentHitContainer->AddHit(layer_id, m_Hit);
        if (m_CurrentShower)
        {
          m_CurrentShower->add_g4hit_id(m_CurrentHitContainer->GetID(), m_Hit->get_hit_id());
        }
        // ownership has been transferred to container, set to null
        if (m_Hit->get_edep() > 0 && (whichactive > 0 || absorbertruth > 0))
        {
          if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
          {
            if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
            {
              pp->SetKeep(1);  // we want to keep the track
            }
          }
        }
        // so we will create a new hit for the next track
        m_Hit = nullptr;
      }
      else
      {
        // if this hit has no energy deposit, just reset it for reuse
        // this means we have to delete it in the dtor. If this was
        // the last hit we processed the memory is still allocated
        m_Hit->Reset();
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4ZDCSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  //now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);
  m_SupportHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_SupportNodeName);
  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4ZDCSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
    gSystem->Exit(1);
  }
  // this is perfectly fine if absorber hits are disabled
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4ZDCSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
    }
  }
  if (!m_SupportHitContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4ZDCSteppingAction::SetTopNode - unable to find " << m_SupportNodeName << std::endl;
    }
  }
}

void PHG4ZDCSteppingAction::SetHitNodeName(const std::string& type, const std::string& name)
{
  if (type == "G4HIT")
  {
    m_HitNodeName = name;
    return;
  }
  else if (type == "G4HIT_ABSORBER")
  {
    m_AbsorberNodeName = name;
    return;
  }
  else if (type == "G4HIT_SUPPORT")
  {
    m_SupportNodeName = name;
    return;
  }
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}

//getting index using copyno
//if in ZDC PMMA fiber
int PHG4ZDCSteppingAction::FindIndexZDC(G4TouchableHandle& touch, int& j, int& k)
{
  G4VPhysicalVolume* envelope = touch->GetVolume(2);  //Get the world
  G4VPhysicalVolume* plate = touch->GetVolume(1);     //Get the fiber plate

  j = envelope->GetCopyNo();
  k = (plate->GetCopyNo()) / 27;

  return 0;
}

int PHG4ZDCSteppingAction::FindIndexSMD(G4TouchableHandle& touch, int& j, int& k)
{
  int jshift = 10;
  int kshift = 10;
  G4VPhysicalVolume* envelope = touch->GetVolume(2);  //Get the envelope
  G4VPhysicalVolume* scint = touch->GetVolume(0);     //Get the fiber plate

  int whichzdc = envelope->GetCopyNo();

  j = scint->GetCopyNo() % 7;
  k = scint->GetCopyNo() / 7;

  if (whichzdc == 1) j += 7;
  // shift the index to avoid conflict with the ZDC index
  j += jshift;
  k += kshift;

  return 0;
}

double PHG4ZDCSteppingAction::ZDCResponse(double beta, double angle)
{
  if (beta < m_BetaThersh) return 0;
  if (angle >= 90) return 0;
  for (int i = 1; i < 9; i++)
  {
    if (beta <= m_Beta[i])
    {
      std::array<double, 18> PMMAsub0 = m_PMMA05[i - 1];
      std::array<double, 18> PMMAsub1 = m_PMMA05[i];
      //find angle bin here and do 1D linear interpolation of angle
      int Abin = (int) angle / 5;
      if (Abin == 0) Abin = 1;
      double avg_ph0 = PMMAsub0[Abin - 1] + (PMMAsub0[Abin] - PMMAsub0[Abin - 1]) * (angle / 5 - Abin + 1);
      double avg_ph1 = PMMAsub1[Abin - 1] + (PMMAsub1[Abin] - PMMAsub1[Abin - 1]) * (angle / 5 - Abin + 1);
      if (avg_ph0 < 0) avg_ph0 = 0;
      if (avg_ph1 < 0) avg_ph1 = 0;
      //linear linear interpolation with beta
      double avg_ph = avg_ph0 + (avg_ph1 - avg_ph0) * (beta - m_Beta[i - 1]) / (m_Beta[i] - m_Beta[i - 1]);
      if (avg_ph < 0) avg_ph = 0;
      //use poisson?
      return avg_ph;
    }
  }

  return 0;
}

double PHG4ZDCSteppingAction::ZDCEResponse(double E, double angle)
{
  if (E < m_E[0]) return 0;

  if (E >= 0.05)
  {
    std::array<double, 36> PMMAsub0 = m_PMMA05E[10];
    int Abin = (int) angle / 5;
    if (Abin == 0) Abin = 1;
    double avg_ph = PMMAsub0[Abin - 1] + (PMMAsub0[Abin] - PMMAsub0[Abin - 1]) * (angle / 5 - Abin + 1);
    return avg_ph;
  }
  else
  {
    for (int i = 1; i < 11; i++)
    {
      if (E <= m_E[i])
      {
        std::array<double, 36> PMMAsub0 = m_PMMA05E[i - 1];
        std::array<double, 36> PMMAsub1 = m_PMMA05E[i];

        int Abin = (int) angle / 5;
        if (Abin == 0) Abin = 1;
        double avg_ph0 = PMMAsub0[Abin - 1] + (PMMAsub0[Abin] - PMMAsub0[Abin - 1]) * (angle / 5 - Abin + 1);
        double avg_ph1 = PMMAsub1[Abin - 1] + (PMMAsub1[Abin] - PMMAsub1[Abin - 1]) * (angle / 5 - Abin + 1);
        if (avg_ph0 < 0) avg_ph0 = 0;
        if (avg_ph1 < 0) avg_ph1 = 0;
        //linear linear interpolation with E
        double avg_ph = avg_ph0 + (avg_ph1 - avg_ph0) * (E - m_E[i - 1]) / (m_E[i] - m_E[i - 1]);
        if (avg_ph < 0) avg_ph = 0;
        //use poisson?
        return avg_ph;
      }
    }
  }
  std::cout << "out of range" << std::endl;
  return 0;
}
