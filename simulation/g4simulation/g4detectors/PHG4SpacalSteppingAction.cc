/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.6 $$
 * \date $$Date: 2015/01/07 23:50:05 $$
 */

#include "PHG4SpacalSteppingAction.h"
#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4SpacalDetector.h"

#include "PHG4CellDefs.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCellGeom_Spacalv1.h"
#include "PHG4CylinderGeom.h"  // for PHG4CylinderGeom
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv1.h"  // for PHG4CylinderGeom_Spaca...

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoDefs.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <phparameter/PHParameters.h>

#include <Geant4/G4IonisParamMat.hh>  // for G4IonisParamMat
#include <Geant4/G4Material.hh>       // for G4Material
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRestD...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>  // for G4TouchableHandle
#include <Geant4/G4Track.hh>            // for G4Track
#include <Geant4/G4TrackStatus.hh>      // for fStopAndKill
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <TSystem.h>

#include <cmath>    // for isfinite
#include <cstdlib>  // for exit
#include <iostream>
#include <string>  // for operator<<, char_traits

class G4VPhysicalVolume;
class PHCompositeNode;

//____________________________________________________________________________..
PHG4SpacalSteppingAction::PHG4SpacalSteppingAction(PHG4SpacalDetector *indetector, const PHParameters *parameters)
  : PHG4SteppingAction(indetector->GetName())
  , m_Detector(indetector)
  , m_Params(parameters)
  , m_doG4Hit(m_Params->get_int_param("saveg4hit"))
  , m_tmin(m_Params->get_double_param("tmin"))
  , m_tmax(m_Params->get_double_param("tmax"))
  , m_dt(m_Params->get_double_param("dt"))
{
}

PHG4SpacalSteppingAction::~PHG4SpacalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

int PHG4SpacalSteppingAction::InitWithNode(PHCompositeNode *topNode)
{
  if (m_doG4Hit)
  {
    return 0;
  }
  PHNodeIterator iter(topNode);
  detector = m_Detector->SuperDetector();
  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // add towerinfo here
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }
  m_CaloInfoContainer = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::EMCAL);
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, "TOWERINFO_SIM_" + detector, "PHObject");
  DetNode->addNode(towerNode);

  return 0;
}

int PHG4SpacalSteppingAction::SetUpGeomNode(PHCompositeNode *topNode)
{
  if (m_geomsetup)
  {
    return 0;
  }

  if (m_doG4Hit)
  {
    return 0;
  }
  PHNodeIterator iter(topNode);
  detector = m_Detector->SuperDetector();

  geonodename = "CYLINDERGEOM_" + detector;
  seggeonodename = "CYLINDERCELLGEOM_" + detector;
  // get the private layergeo
  _layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename);
  if (!_layergeo)
  {
    std::cout << "PHG4SpacalSteppingAction::InitWithNode - Fatal Error - Could not locate sim geometry node "
              << geonodename << std::endl;
    exit(1);
  }

  _seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename);
  if (!_seggeo)
  {
    std::cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - could not locate cell geometry node "
              << seggeonodename << std::endl;
    exit(1);
  }
  PHG4CylinderCellGeom *_geo_raw = _seggeo->GetFirstLayerCellGeom();
  _geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(_geo_raw);
  assert(_geo);
  const PHG4CylinderGeom *_layergeom_raw = _layergeo->GetFirstLayerGeom();
  assert(_layergeom_raw);
  // a special implimentation of PHG4CylinderGeom is required here.
  _layergeom = dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(_layergeom_raw);
  assert(_layergeom);
  m_geomsetup = true;
  return 0;
}

bool PHG4SpacalSteppingAction::NoHitSteppingAction(const G4Step *aStep)
{
  // get volume of the current step
  G4VPhysicalVolume *volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  int isactive = m_Detector->IsInCylinderActive(volume);
  if (isactive > PHG4SpacalDetector::INACTIVE)
  {
    G4StepPoint *prePoint = aStep->GetPreStepPoint();
    G4StepPoint *postPoint = aStep->GetPostStepPoint();
    // time window cut
    double pretime = prePoint->GetGlobalTime() / nanosecond;
    double posttime = postPoint->GetGlobalTime() / nanosecond;
    if (posttime < m_tmin || pretime > m_tmax)
    {
      return false;
    }
    if ((posttime - pretime) > m_dt)
    {
      return false;
    }

    int scint_id = -1;

    if (  //
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_SameLengthFiberPerTower ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_Tilted ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower  //
    )
    {
      // SPACAL ID that is associated with towers
      int sector_ID = 0;
      int tower_ID = 0;
      int fiber_ID = 0;

      if (isactive == PHG4SpacalDetector::FIBER_CORE)
      {
        fiber_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
        tower_ID = prePoint->GetTouchable()->GetReplicaNumber(2);
        sector_ID = prePoint->GetTouchable()->GetReplicaNumber(3);
      }

      else
      {
        return false;
      }

      // compact the tower/sector/fiber ID into 32 bit scint_id, so we could save some space for SPACAL hits
      scint_id = PHG4CylinderGeom_Spacalv3::scint_id_coder(sector_ID, tower_ID, fiber_ID).scint_ID;
    }
    else
    {
      // other configuraitons
      if (isactive == PHG4SpacalDetector::FIBER_CORE)
      {
        scint_id = prePoint->GetTouchable()->GetReplicaNumber(2);
      }
      else
      {
        return false;
      }
    }
    // decode scint_id
    PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);

    // convert to z_ID, phi_ID
    std::pair<int, int> tower_z_phi_ID = _layergeom->get_tower_z_phi_ID(decoder.tower_ID, decoder.sector_ID);
    const int &tower_ID_z = tower_z_phi_ID.first;
    const int &tower_ID_phi = tower_z_phi_ID.second;

    PHG4CylinderGeom_Spacalv3::tower_map_t::const_iterator it_tower =
        _layergeom->get_sector_tower_map().find(decoder.tower_ID);
    assert(it_tower != _layergeom->get_sector_tower_map().end());

    // convert tower_ID_z to to eta bin number
    int etabin = -1;
    try
    {
      etabin = _geo->get_etabin_block(tower_ID_z);  // block eta bin
    }
    catch (std::exception &e)
    {
      std::cout << "Print cell geometry:" << std::endl;
      _geo->identify();
      std::cout << "Print scint_id_coder:" << std::endl;
      decoder.identify();
      std::cout << "PHG4SpacalSteppingAction::UserSteppingAction::"
                << " - Fatal Error - " << e.what() << std::endl;
      exit(1);
    }

    const int sub_tower_ID_x = it_tower->second.get_sub_tower_ID_x(decoder.fiber_ID);
    const int sub_tower_ID_y = it_tower->second.get_sub_tower_ID_y(decoder.fiber_ID);
    unsigned short etabinshort = etabin * _layergeom->get_n_subtower_eta() + sub_tower_ID_y;
    unsigned short phibin = tower_ID_phi * _layergeom->get_n_subtower_phi() + sub_tower_ID_x;

    // get light yield
    double light_yield = GetVisibleEnergyDeposition(aStep);

    if (light_collection_model.use_fiber_model())
    {
      const G4TouchableHandle &theTouchable0 = prePoint->GetTouchableHandle();
      const G4ThreeVector &worldPosition0 = prePoint->GetPosition();
      G4ThreeVector localPosition = theTouchable0->GetHistory()->GetTopTransform().TransformPoint(worldPosition0);
      const double localz0 = localPosition.z();
      // post point
      const G4TouchableHandle &theTouchable1 = postPoint->GetTouchableHandle();
      const G4ThreeVector &worldPosition1 = postPoint->GetPosition();
      localPosition = theTouchable1->GetHistory()->GetTopTransform().TransformPoint(worldPosition1);
      const double localz1 = localPosition.z();

      const double z = 0.5 * (localz0 + localz1);
      assert(not std::isnan(z));

      light_yield *= light_collection_model.get_fiber_transmission(z);
    }

    // light yield correction from light guide collection efficiency:
    if (light_collection_model.use_fiber_model())
    {
      const double x = it_tower->second.get_position_fraction_x_in_sub_tower(decoder.fiber_ID);
      const double y = it_tower->second.get_position_fraction_y_in_sub_tower(decoder.fiber_ID);

      light_yield *= light_collection_model.get_light_guide_efficiency(x, y);
    }
    unsigned int tower_key = TowerInfoDefs::encode_emcal(etabinshort, phibin);
    m_CaloInfoContainer->get_tower_at_key(tower_key)->set_energy(m_CaloInfoContainer->get_tower_at_key(tower_key)->get_energy() + light_yield);

    // set keep for the track
    const G4Track *aTrack = aStep->GetTrack();
    if (light_yield > 0)
    {
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
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
bool PHG4SpacalSteppingAction::UserSteppingAction(const G4Step *aStep, bool)
{
  if (!m_doG4Hit)
  {
    return NoHitSteppingAction(aStep);
  }

  // get volume of the current step
  G4VPhysicalVolume *volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  const G4Track *aTrack = aStep->GetTrack();

  int layer_id = m_Detector->get_Layer();
  // make sure we are in a volume
  // IsInCylinderActive returns the number of the scintillator
  // slat which has fired
  int isactive = m_Detector->IsInCylinderActive(volume);
  if (isactive > PHG4SpacalDetector::INACTIVE)
  {
    bool geantino = false;
    // the check for the pdg code speeds things up, I do not want to make
    // an expensive string compare for every track when we know
    // geantino or chargedgeantino has pid=0
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
    {
      geantino = true;
    }
    G4StepPoint *prePoint = aStep->GetPreStepPoint();
    G4StepPoint *postPoint = aStep->GetPostStepPoint();
    int scint_id = -1;

    if (  //
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_SameLengthFiberPerTower ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_Tilted ||
        m_Detector->get_geom()->get_config() == PHG4SpacalDetector::SpacalGeom_t::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower  //
    )
    {
      // SPACAL ID that is associated with towers
      int sector_ID = 0;
      int tower_ID = 0;
      int fiber_ID = 0;

      if (isactive == PHG4SpacalDetector::FIBER_CORE)
      {
        fiber_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
        tower_ID = prePoint->GetTouchable()->GetReplicaNumber(2);
        sector_ID = prePoint->GetTouchable()->GetReplicaNumber(3);
      }

      else if (isactive == PHG4SpacalDetector::FIBER_CLADING)
      {
        fiber_ID = prePoint->GetTouchable()->GetReplicaNumber(0);
        tower_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
        sector_ID = prePoint->GetTouchable()->GetReplicaNumber(2);
      }

      else if (isactive == PHG4SpacalDetector::ABSORBER)
      {
        tower_ID = prePoint->GetTouchable()->GetReplicaNumber(0);
        sector_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
      }

      else if (isactive == PHG4SpacalDetector::SUPPORT)
      {
        tower_ID = prePoint->GetTouchable()->GetReplicaNumber(0);
        sector_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
        fiber_ID = (1 << (PHG4CylinderGeom_Spacalv3::scint_id_coder::kfiber_bit)) - 1;  // use max fiber ID to flag for support strucrtures.

        //        std::cout <<"PHG4SpacalSteppingAction::UserSteppingAction - SUPPORT tower_ID = "<<tower_ID<<std::endl;
      }

      // compact the tower/sector/fiber ID into 32 bit scint_id, so we could save some space for SPACAL hits
      scint_id = PHG4CylinderGeom_Spacalv3::scint_id_coder(sector_ID, tower_ID, fiber_ID).scint_ID;
    }
    else
    {
      // other configuraitons
      if (isactive == PHG4SpacalDetector::FIBER_CORE)
      {
        scint_id = prePoint->GetTouchable()->GetReplicaNumber(2);
      }
      else if (isactive == PHG4SpacalDetector::FIBER_CLADING)
      {
        scint_id = prePoint->GetTouchable()->GetReplicaNumber(1);
      }
      else
      {
        scint_id = prePoint->GetTouchable()->GetReplicaNumber(0);
      }
    }

    //       std::cout << "track id " << aTrack->GetTrackID() << std::endl;
    //        std::cout << "time prepoint: " << prePoint->GetGlobalTime() << std::endl;
    //        std::cout << "time postpoint: " << postPoint->GetGlobalTime() << std::endl;
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      // if previous hit was saved, hit pointer was set to nullptr
      // and we have to make a new one
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }
      m_Hit->set_layer((unsigned int) layer_id);
      m_Hit->set_scint_id(scint_id);  // isactive contains the scintillator slat id
      // here we set the entrance values in cm
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      // time in ns
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      // set the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      m_SaveTrackid = aTrack->GetTrackID();
      // set the initial energy deposit
      m_Hit->set_edep(0);
      // Now add the hit
      if (isactive == PHG4SpacalDetector::FIBER_CORE)  // the slat ids start with zero
      {
        // store all pre local coordinates
        StoreLocalCoordinate(m_Hit, aStep, true, false);
        m_Hit->set_eion(0);  // only implemented for v5 otherwise empty
        m_Hit->set_light_yield(0);
        m_CurrentHitContainer = m_HitContainer;
      }
      else
      {
        m_CurrentHitContainer = m_AbsorberHitContainer;
      }
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_CurrentShower = pp->GetShower();
        }
      }

      if (m_Hit->get_z(0) > get_zmax() || m_Hit->get_z(0) < get_zmin())
      {
        std::cout << "PHG4SpacalSteppingAction: hit outside acceptance, layer: "
                  << layer_id << std::endl;
        m_Hit->identify();
      }
      break;
    default:
      break;
    }
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!m_Hit || !std::isfinite(m_Hit->get_x(0)))
    {
      std::cout << GetName() << ": hit was not created" << std::endl;
      std::cout << "prestep status: " << prePoint->GetStepStatus()
                << ", last post step status: " << m_SavePostStepStatus << std::endl;
      exit(1);
    }
    m_SavePostStepStatus = postPoint->GetStepStatus();
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != m_SaveTrackid)
    {
      std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
      std::cout << "saved track: " << m_SaveTrackid
                << ", current trackid: " << aTrack->GetTrackID()
                << std::endl;
      exit(1);
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
    // sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);

    if (isactive == PHG4SpacalDetector::FIBER_CORE)  // only for active areas
    {
      // store all pre local coordinates
      StoreLocalCoordinate(m_Hit, aStep, false, true);

      m_Hit->set_eion(m_Hit->get_eion() + eion);

      double light_yield = GetVisibleEnergyDeposition(aStep);

      static bool once = true;
      if (once and edep > 0)
      {
        once = false;

        if (Verbosity() > 0)
        {
          std::cout << "PHG4SpacalSteppingAction::UserSteppingAction::"
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

      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    }

    if (m_Hit->get_z(1) > get_zmax() || m_Hit->get_z(1) < get_zmin())
    {
      std::cout << "PHG4SpacalSteppingAction: hit outside acceptance get_zmin() "
                << get_zmin() << ", get_zmax() " << get_zmax() << " at exit"
                << std::endl;
      m_Hit->identify();
    }
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
                            //          m_Hit->set_eion(-1);
    }
    if (edep > 0)
    {
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
    }
    // if any of these conditions is true this is the last step in
    // this volume and we need to save the hit
    // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
    // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
    // (happens when your detector goes outside world volume)
    // aTrack->GetTrackStatus() == fStopAndKill: track ends
    if (postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (m_Hit->get_edep())
      {
        m_CurrentHitContainer->AddHit(layer_id, m_Hit);
        if (m_CurrentShower)
        {
          m_CurrentShower->add_g4hit_id(m_CurrentHitContainer->GetID(), m_Hit->get_hit_id());
        }
        // ownership has been transferred to container, set to null
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
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4SpacalSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  // we can only play with the geometry node after ConstructMe is called
  if ((!m_geomsetup) && (!m_doG4Hit))
  {
    SetUpGeomNode(topNode);
  }
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);
  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4SpacalSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
    gSystem->Exit(1);
  }
  // this is perfectly fine if absorber hits are disabled
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4SpacalSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
    }
  }
}

double
PHG4SpacalSteppingAction::get_zmin() const
{
  if (!m_Detector)
  {
    return 0;
  }
  else
  {
    return m_Detector->get_geom()->get_zmin() - .0001;
  }
}

double
PHG4SpacalSteppingAction::get_zmax() const
{
  if (!m_Detector)
  {
    return 0;
  }
  else
  {
    return m_Detector->get_geom()->get_zmax() + .0001;
  }
}

void PHG4SpacalSteppingAction::SetHitNodeName(const std::string &type, const std::string &name)
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
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}
