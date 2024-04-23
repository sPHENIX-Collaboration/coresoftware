#include "PHG4Reco.h"

#include "Fun4AllMessenger.h"
#include "G4TBMagneticFieldSetup.hh"
#include "PHG4DisplayAction.h"
#include "PHG4InEvent.h"
#include "PHG4PhenixDetector.h"
#include "PHG4PhenixDisplayAction.h"
#include "PHG4PhenixEventAction.h"
#include "PHG4PhenixStackingAction.h"
#include "PHG4PhenixSteppingAction.h"
#include "PHG4PhenixTrackingAction.h"
#include "PHG4PrimaryGeneratorAction.h"
#include "PHG4Subsystem.h"
#include "PHG4TrackingAction.h"
#include "PHG4UIsession.h"
#include "PHG4Utils.h"

#include <g4decayer/EDecayType.hh>
#include <g4decayer/P6DExtDecayerPhysics.hh>

#include <g4decayer/EvtGenExtDecayerPhysics.hh>

#include <phgeom/PHGeomUtility.h>

#include <g4gdml/PHG4GDMLUtility.hh>

#include <phfield/PHFieldConfig.h>  // for PHFieldConfig
#include <phfield/PHFieldConfigv1.h>
#include <phfield/PHFieldConfigv2.h>
#include <phfield/PHFieldUtility.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <TSystem.h>  // for TSystem, gSystem

#include <CLHEP/Random/Random.h>

#include <G4HadronicParameters.hh>  // for G4HadronicParameters
#include <Geant4/G4Cerenkov.hh>
#include <Geant4/G4Element.hh>       // for G4Element
#include <Geant4/G4EventManager.hh>  // for G4EventManager
#include <Geant4/G4HadronicProcessStore.hh>
#include <Geant4/G4IonisParamMat.hh>  // for G4IonisParamMat
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4OpAbsorption.hh>
#include <Geant4/G4OpBoundaryProcess.hh>
#include <Geant4/G4OpMieHG.hh>
#include <Geant4/G4OpRayleigh.hh>
#include <Geant4/G4OpWLS.hh>
#include <Geant4/G4OpticalPhoton.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4PhotoElectricEffect.hh>  // for G4PhotoElectricEffect
#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4RunManager.hh>
#include <Geant4/G4Scintillation.hh>
#include <Geant4/G4StepLimiterPhysics.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>  // for G4double, G4int
#include <Geant4/G4UIExecutive.hh>
#include <Geant4/G4UImanager.hh>
#include <Geant4/G4UImessenger.hh>          // for G4UImessenger
#include <Geant4/G4VModularPhysicsList.hh>  // for G4VModularPhysicsList
#include <Geant4/G4Version.hh>
#include <Geant4/G4VisExecutive.hh>
#include <Geant4/G4VisManager.hh>  // for G4VisManager
#include <Geant4/Randomize.hh>     // for G4Random

// physics lists
#include <Geant4/FTFP_BERT.hh>
#include <Geant4/FTFP_BERT_HP.hh>
#include <Geant4/FTFP_INCLXX.hh>
#include <Geant4/FTFP_INCLXX_HP.hh>
#include <Geant4/QGSP_BERT.hh>
#include <Geant4/QGSP_BERT_HP.hh>
#include <Geant4/QGSP_BIC.hh>
#include <Geant4/QGSP_BIC_HP.hh>
#include <Geant4/QGSP_INCLXX.hh>
#include <Geant4/QGSP_INCLXX_HP.hh>

#include <cassert>
#include <cstdlib>
#include <exception>  // for exception
#include <filesystem>
#include <iostream>   // for operator<<, endl
#include <memory>

class G4EmSaturation;
class G4TrackingManager;
class G4VPhysicalVolume;
class PHField;
class PHG4EventAction;
class PHG4StackingAction;
class PHG4SteppingAction;

//_________________________________________________________________
PHG4Reco::PHG4Reco(const std::string &name)
  : SubsysReco(name)
  , m_Fun4AllMessenger(new Fun4AllMessenger(Fun4AllServer::instance()))
{
  return;
}

//_________________________________________________________________
PHG4Reco::~PHG4Reco()
{
  // one can delete null pointer (it results in a nop), so checking if
  // they are non zero is not needed
  delete m_Field;
  delete m_RunManager;
  delete m_UISession;
  delete m_VisManager;
  delete m_Fun4AllMessenger;
  while (m_SubsystemList.begin() != m_SubsystemList.end())
  {
    delete m_SubsystemList.back();
    m_SubsystemList.pop_back();
  }
  delete m_DisplayAction;
}

//_________________________________________________________________
int PHG4Reco::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "========================= PHG4Reco::Init() ================================" << std::endl;
  }
  unsigned int iseed = PHRandomSeed();
  G4Seed(iseed);  // fixed seed handled in PHRandomSeed()

  // create GEANT run manager
  if (Verbosity() > 1)
  {
    std::cout << "PHG4Reco::Init - create run manager" << std::endl;
  }

  // redirect GEANT verbosity to nowhere
  //  if (Verbosity() < 1)
  if (false)
  {
    G4UImanager *uimanager = G4UImanager::GetUIpointer();
    m_UISession = new PHG4UIsession();
    uimanager->SetSession(m_UISession);
    uimanager->SetCoutDestination(m_UISession);
  }

  m_RunManager = new G4RunManager();

  DefineMaterials();
  // create physics processes
  G4VModularPhysicsList *myphysicslist = nullptr;
  if (m_PhysicsList == "FTFP_BERT")
  {
    myphysicslist = new FTFP_BERT(Verbosity());
  }
  else if (m_PhysicsList == "FTFP_BERT_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new FTFP_BERT_HP(Verbosity());
  }
  else if (m_PhysicsList == "FTFP_INCLXX")
  {
    myphysicslist = new FTFP_INCLXX(Verbosity());
  }
  else if (m_PhysicsList == "FTFP_INCLXX_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new FTFP_INCLXX_HP(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_BERT")
  {
    myphysicslist = new QGSP_BERT(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_BERT_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new QGSP_BERT_HP(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_BIC")
  {
    myphysicslist = new QGSP_BIC(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_BIC_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new QGSP_BIC_HP(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_INCLXX")
  {
    myphysicslist = new QGSP_INCLXX(Verbosity());
  }
  else if (m_PhysicsList == "QGSP_INCLXX_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new QGSP_INCLXX_HP(Verbosity());
  }
  else
  {
    std::cout << "Physics List " << m_PhysicsList << " not implemented" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  if (m_Decayer == kPYTHIA6Decayer)
  {
    std::cout << "Use PYTHIA Decayer" << std::endl;
    G4HadronicParameters::Instance()->SetEnableBCParticles(false);  // Disable the Geant4 built in HF Decay and use external decayers for them
    P6DExtDecayerPhysics *decayer = new P6DExtDecayerPhysics();
    if (m_ActiveForceDecayFlag)
    {
      decayer->SetForceDecay(m_ForceDecayType);
    }
    myphysicslist->RegisterPhysics(decayer);
  }

  if (m_Decayer == kEvtGenDecayer)
  {
    std::cout << "Use EvtGen Decayer" << std::endl;
    G4HadronicParameters::Instance()->SetEnableBCParticles(false);  // Disable the Geant4 built in HF Decay and use external decayers for them
    EvtGenExtDecayerPhysics *decayer = new EvtGenExtDecayerPhysics();
    if (CustomizeDecay)
    {
      decayer->CustomizeEvtGenDecay(EvtGenDecayFile);
    }

    myphysicslist->RegisterPhysics(decayer);
  }

  if (m_Decayer == kGEANTInternalDecayer)
  {
    std::cout << "Use GEANT Internal Decayer" << std::endl;
  }

  myphysicslist->RegisterPhysics(new G4StepLimiterPhysics());
  // initialize cuts so we can ask the world region for it's default
  // cuts to propagate them to other regions in DefineRegions()
  myphysicslist->SetCutsWithDefault();
  m_RunManager->SetUserInitialization(myphysicslist);

  DefineRegions();
  // initialize registered subsystems
  for (SubsysReco *reco : m_SubsystemList)
  {
    reco->Init(topNode);
  }

  if (Verbosity() > 0)
  {
    std::cout << "===========================================================================" << std::endl;
  }
  return 0;
}

int PHG4Reco::InitField(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHG4Reco::InitField - create magnetic field setup" << std::endl;
  }

  std::unique_ptr<PHFieldConfig> default_field_cfg(nullptr);

  if (std::filesystem::path(m_FieldMapFile).extension() != ".root")
  {
    // loading from database
    std::string url = CDBInterface::instance()->getUrl(m_FieldMapFile);
    default_field_cfg.reset(new PHFieldConfigv1(m_FieldConfigType, url, m_MagneticFieldRescale));
  }
  else if (m_FieldMapFile != "NONE")
  {
    default_field_cfg.reset(new PHFieldConfigv1(m_FieldConfigType, m_FieldMapFile, m_MagneticFieldRescale));
  }
  else
  {
    default_field_cfg.reset(new PHFieldConfigv2(0, 0, m_MagneticField * m_MagneticFieldRescale));
  }

  if (Verbosity() > 1)
  {
    std::cout << "PHG4Reco::InitField - create magnetic field setup" << std::endl;
  }

  PHField *phfield = PHFieldUtility::GetFieldMapNode(default_field_cfg.get(), topNode, Verbosity() + 1);
  assert(phfield);

  m_Field = new G4TBMagneticFieldSetup(phfield);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4Reco::InitRun(PHCompositeNode *topNode)
{
  // this is a dumb protection against executing this twice.
  // we have cases (currently detector display or material scan) where
  // we need the detector but have not run any events (who wants to wait
  // for processing an event if you just want a detector display?).
  // Then the InitRun is executed from a macro. But if you decide to run events
  // afterwards, the InitRun is executed by the framework together with all
  // other modules. This prevents the PHG4Reco::InitRun() to be executed
  // again in this case
  static int InitRunDone = 0;
  if (InitRunDone)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  InitRunDone = 1;
  if (Verbosity() > 0)
  {
    std::cout << "========================= PHG4Reco::InitRun() ================================" << std::endl;
  }

  recoConsts *rc = recoConsts::instance();

  rc->set_StringFlag("WorldMaterial", m_WorldMaterial);
  // build world material - so in subsequent code we can call
  //  G4Material::GetMaterial(rc->get_StringFlag("WorldMaterial"))
  // if the world material is not in the nist DB, we need to implement it here
  G4NistManager::Instance()->FindOrBuildMaterial(m_WorldMaterial);
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_Be");
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  // G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  rc->set_FloatFlag("WorldSizex", m_WorldSize[0]);
  rc->set_FloatFlag("WorldSizey", m_WorldSize[1]);
  rc->set_FloatFlag("WorldSizez", m_WorldSize[2]);

  // setup the global field
  const int field_ret = InitField(topNode);
  if (field_ret != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << "PHG4Reco::InitRun- Error - Failed field init with status = " << field_ret << std::endl;
    return field_ret;
  }

  // initialize registered subsystems
  for (SubsysReco *reco : m_SubsystemList)
  {
    if (Verbosity() >= 1)
    {
      std::cout << "PHG4Reco::InitRun - " << reco->Name() << "->InitRun" << std::endl;
    }
    reco->InitRun(topNode);
  }

  // create phenix detector, add subsystems, and register to GEANT
  // create display settings before detector
  m_DisplayAction = new PHG4PhenixDisplayAction(Name());
  if (Verbosity() > 1)
  {
    std::cout << "PHG4Reco::Init - create detector" << std::endl;
  }
  m_Detector = new PHG4PhenixDetector(this);
  m_Detector->Verbosity(Verbosity());
  m_Detector->SetWorldSizeX(m_WorldSize[0] * cm);
  m_Detector->SetWorldSizeY(m_WorldSize[1] * cm);
  m_Detector->SetWorldSizeZ(m_WorldSize[2] * cm);
  m_Detector->SetWorldShape(m_WorldShape);
  m_Detector->SetWorldMaterial(m_WorldMaterial);

  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    if (g4sub->GetDetector())
    {
      m_Detector->AddDetector(g4sub->GetDetector());
    }
  }
  m_RunManager->SetUserInitialization(m_Detector);

  if (m_disableUserActions)
  {
    std::cout << "PHG4Reco::InitRun - WARNING - event/track/stepping action disabled! "
              << "This is aimed to reduce resource consumption for G4 running only. E.g. dose analysis. "
              << "Meanwhile, it will disable all Geant4 based analysis. Toggle this feature on/off with PHG4Reco::setDisableUserActions()" << std::endl;
  }

  setupInputEventNodeReader(topNode);
  // create main event action, add subsystemts and register to GEANT
  m_EventAction = new PHG4PhenixEventAction();

  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    PHG4EventAction *evtact = g4sub->GetEventAction();
    if (evtact)
    {
      m_EventAction->AddAction(evtact);
    }
  }

  if (not m_disableUserActions)
  {
    m_RunManager->SetUserAction(m_EventAction);
  }

  // create main stepping action, add subsystems and register to GEANT
  m_StackingAction = new PHG4PhenixStackingAction();
  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    PHG4StackingAction *action = g4sub->GetStackingAction();
    if (action)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Adding steppingaction for " << g4sub->Name() << std::endl;
      }
      m_StackingAction->AddAction(g4sub->GetStackingAction());
    }
  }

  if (not m_disableUserActions)
  {
    m_RunManager->SetUserAction(m_StackingAction);
  }

  // create main stepping action, add subsystems and register to GEANT
  m_SteppingAction = new PHG4PhenixSteppingAction();
  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    PHG4SteppingAction *action = g4sub->GetSteppingAction();
    if (action)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Adding steppingaction for " << g4sub->Name() << std::endl;
      }

      m_SteppingAction->AddAction(g4sub->GetSteppingAction());
    }
  }

  if (not m_disableUserActions)
  {
    m_RunManager->SetUserAction(m_SteppingAction);
  }

  // create main tracking action, add subsystems and register to GEANT
  m_TrackingAction = new PHG4PhenixTrackingAction();
  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    m_TrackingAction->AddAction(g4sub->GetTrackingAction());

    // not all subsystems define a user tracking action
    if (g4sub->GetTrackingAction())
    {
      // make tracking manager accessible within user tracking action if defined
      if (G4TrackingManager *trackingManager = G4EventManager::GetEventManager()->GetTrackingManager())
      {
        g4sub->GetTrackingAction()->SetTrackingManagerPointer(trackingManager);
      }
    }
  }

  if (not m_disableUserActions)
  {
    m_RunManager->SetUserAction(m_TrackingAction);
  }

  // initialize
  m_RunManager->Initialize();

#if G4VERSION_NUMBER >= 1033
  G4EmSaturation *emSaturation = G4LossTableManager::Instance()->EmSaturation();
  if (!emSaturation)
  {
    std::cout << PHWHERE << "Could not initialize EmSaturation, Birks constants will fail" << std::endl;
  }
#endif

  // add cerenkov and optical photon processes
  // std::cout << std::endl << "Ignore the next message - we implemented this correctly" << std::endl;
  G4Cerenkov *theCerenkovProcess = new G4Cerenkov("Cerenkov");
  // std::cout << "End of bogus warning message" << std::endl << std::endl;
  G4Scintillation *theScintillationProcess = new G4Scintillation("Scintillation");

  /*
    if (Verbosity() > 0)
    {
    // This segfaults
    theCerenkovProcess->DumpPhysicsTable();
    }
  */
  theCerenkovProcess->SetMaxNumPhotonsPerStep(300);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(false);  // current PHG4TruthTrackingAction does not support suspect active track and track secondary first
#if G4VERSION_NUMBER < 1100
  theScintillationProcess->SetScintillationYieldFactor(1.0);
#endif
  theScintillationProcess->SetTrackSecondariesFirst(false);
  // theScintillationProcess->SetScintillationExcitationRatio(1.0);

  // Use Birks Correction in the Scintillation process

  // G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  // theScintillationProcess->AddSaturation(emSaturation);

  G4ParticleTable *theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator *_theParticleIterator;
  _theParticleIterator = theParticleTable->GetIterator();
  _theParticleIterator->reset();
  while ((*_theParticleIterator)())
  {
    G4ParticleDefinition *particle = _theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (theCerenkovProcess->IsApplicable(*particle))
    {
      pmanager->AddProcess(theCerenkovProcess);
      pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
    }
    if (theScintillationProcess->IsApplicable(*particle))
    {
      pmanager->AddProcess(theScintillationProcess);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    }
    for (PHG4Subsystem *g4sub : m_SubsystemList)
    {
      g4sub->AddProcesses(particle);
    }
  }
  G4ProcessManager *pmanager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  // std::cout << " AddDiscreteProcess to OpticalPhoton " << std::endl;
  pmanager->AddDiscreteProcess(new G4OpAbsorption());
  pmanager->AddDiscreteProcess(new G4OpRayleigh());
  pmanager->AddDiscreteProcess(new G4OpMieHG());
  pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
  pmanager->AddDiscreteProcess(new G4OpWLS());
  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
  // pmanager->DumpInfo();

  // needs large amount of memory which kills central hijing events
  // store generated trajectories
  // if( G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager() ){
  //  trackingManager->SetStoreTrajectory( true );
  //}

  // quiet some G4 print-outs (EM and Hadronic settings during first event)
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4LossTableManager::Instance()->SetVerbose(1);

  if ((Verbosity() < 1) && (m_UISession))
  {
    m_UISession->Verbosity(1);  // let messages after setup come through
  }

  // Geometry export to DST
  if (m_SaveDstGeometryFlag)
  {
    const std::string filename = PHGeomUtility::GenerateGeometryFileName("gdml");
    std::cout << "PHG4Reco::InitRun - export geometry to DST via tmp file " << filename << std::endl;

    Dump_GDML(filename);

    PHGeomUtility::ImportGeomFile(topNode, filename);

    PHGeomUtility::RemoveGeometryFile(filename);
  }

  if (Verbosity() > 0)
  {
    std::cout << "===========================================================================" << std::endl;
  }

  // dump geometry to root file
  if (m_ExportGeometry)
  {
    std::cout << "PHG4Reco::InitRun - writing geometry to " << m_ExportGeomFilename << std::endl;
    PHGeomUtility::ExportGeomtry(topNode, m_ExportGeomFilename);
  }

  if (PHRandomSeed::Verbosity() >= 2)
  {
    // at high verbosity, to save the random number to file
    G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  }
  return 0;
}

//________________________________________________________________
// Dump TGeo File
void PHG4Reco::Dump_GDML(const std::string &filename)
{
  PHG4GDMLUtility ::Dump_GDML(filename, m_Detector->GetPhysicalVolume());
}

//________________________________________________________________
// Dump TGeo File using native Geant4 tools
void PHG4Reco::Dump_G4_GDML(const std::string &filename)
{
  PHG4GDMLUtility::Dump_G4_GDML(filename, m_Detector->GetPhysicalVolume());
}

//_________________________________________________________________
int PHG4Reco::ApplyCommand(const std::string &cmd)
{
  InitUImanager();
  int iret = m_UImanager->ApplyCommand(cmd.c_str());
  return iret;
}
//_________________________________________________________________

int PHG4Reco::StartGui()
{
  // kludge, using boost::dll::program_location().string().c_str() for the
  // program name and putting it into args lead to invalid reads in G4String
  char *args[] = {(char *) ("root.exe")};
  G4UIExecutive *ui = new G4UIExecutive(1, args, "qt");
  InitUImanager();
  m_UImanager->ApplyCommand("/control/execute init_gui_vis.mac");
  ui->SessionStart();
  delete ui;
  return 0;
}

int PHG4Reco::InitUImanager()
{
  if (!m_UImanager)
  {
    // Get the pointer to the User Interface manager
    // Initialize visualization
    m_VisManager = new G4VisExecutive;
    m_VisManager->Initialize();
    m_UImanager = G4UImanager::GetUIpointer();
  }
  return 0;
}

//_________________________________________________________________
int PHG4Reco::process_event(PHCompositeNode *topNode)
{
  if (PHRandomSeed::Verbosity() >= 2)
  {
    G4Random::showEngineStatus();
  }

  // make sure Actions and subsystems have the relevant pointers set
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  m_GeneratorAction->SetInEvent(ineve);

  for (SubsysReco *reco : m_SubsystemList)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "PHG4Reco::process_event - " << reco->Name() << "->process_event" << std::endl;
    }

    try
    {
      reco->process_event(topNode);
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " caught exception thrown during process_event from "
                << reco->Name() << std::endl;
      std::cout << "error: " << e.what() << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // run one event
  if (Verbosity() >= 2)
  {
    std::cout << " PHG4Reco::process_event - "
              << "run one event :" << std::endl;
    ineve->identify();
  }
  m_RunManager->BeamOn(1);

  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    if (Verbosity() >= 2)
    {
      std::cout << " PHG4Reco::process_event - " << g4sub->Name() << "->process_after_geant" << std::endl;
    }
    try
    {
      g4sub->process_after_geant(topNode);
    }
    catch (const std::exception &e)
    {
      std::cout << PHWHERE << " caught exception thrown during process_after_geant from "
                << g4sub->Name() << std::endl;
      std::cout << "error: " << e.what() << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  return 0;
}

int PHG4Reco::ResetEvent(PHCompositeNode *topNode)
{
  for (SubsysReco *reco : m_SubsystemList)
  {
    reco->ResetEvent(topNode);
  }
  return 0;
}

void PHG4Reco::Print(const std::string &what) const
{
  for (SubsysReco *reco : m_SubsystemList)
  {
    if (what.empty() || what == "ALL" || (reco->Name()).find(what) != std::string::npos)
    {
      std::cout << "Printing " << reco->Name() << std::endl;
      reco->Print(what);
      std::cout << "---------------------------" << std::endl;
    }
  }
  return;
}

int PHG4Reco::setupInputEventNodeReader(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    ineve = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }
  // check if we have already registered a generator before creating the default which uses PHG4InEvent Node
  if (!m_GeneratorAction)
  {
    m_GeneratorAction = new PHG4PrimaryGeneratorAction();
  }
  m_RunManager->SetUserAction(m_GeneratorAction);
  return 0;
}

void PHG4Reco::set_rapidity_coverage(const double eta)
{
  m_EtaCoverage = eta;
  PHG4Utils::SetPseudoRapidityCoverage(eta);
}

void PHG4Reco::G4Seed(const unsigned int i)
{
  CLHEP::HepRandom::setTheSeed(i);

  if (PHRandomSeed::Verbosity() >= 2)
  {
    G4Random::showEngineStatus();
  }

  return;
}

//____________________________________________________________________________
void PHG4Reco::DefineMaterials()
{
  G4String symbol, name;  // a=mass of a mole;
  G4double density;       // z=mean number of protons;
  G4double fractionmass, a;
  G4int ncomponents, natoms, z;
  // this is for FTFP_BERT_HP where the neutron code barfs
  // if the z difference to the last known element (U) is too large
  // home made compounds
  // this is a legacy implementation but if they are used in multiple
  // subsystems put them here
  // otherwise implement locally used ones now in the DefineMaterials()
  // method in your subsystem

  // making quartz
  G4Material *quartz = new G4Material("Quartz", density = 2.200 * g / cm3, ncomponents = 2);
  quartz->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 1);
  quartz->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 2);

  // making carbon fiber epoxy
  G4Material *cfrp_intt = new G4Material("CFRP_INTT", density = 1.69 * g / cm3, ncomponents = 3);
  cfrp_intt->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 10);
  cfrp_intt->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 6);
  cfrp_intt->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 1);

  // water glycol mixture for the INTT endcap rings
  G4Material *PropyleneGlycol = new G4Material("Propyleneglycol", 1.036 * g / cm3, 3);
  PropyleneGlycol->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 3);
  PropyleneGlycol->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 8);
  PropyleneGlycol->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 2);

  G4Material *WaterGlycol_INTT = new G4Material("WaterGlycol_INTT", density = (0.997 * 0.7 + 1.036 * 0.3) * g / cm3, ncomponents = 2);
  WaterGlycol_INTT->AddMaterial(PropyleneGlycol, fractionmass = 0.30811936);
  WaterGlycol_INTT->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER"), fractionmass = 0.69188064);

  // making Rohacell foam 110
  G4Material *rohacell_foam_110 = new G4Material("ROHACELL_FOAM_110", density = 0.110 * g / cm3, ncomponents = 4);
  rohacell_foam_110->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 8);
  rohacell_foam_110->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 11);
  rohacell_foam_110->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 2);
  rohacell_foam_110->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), 1);

  // making Rohacell foam ROHACELL 51 WF
  // Source of density: https://www.rohacell.com/product/peek-industrial/downloads/rohacell%20wf%20product%20information.pdf
  G4Material *rohacell_foam_51 = new G4Material("ROHACELL_FOAM_51", density = 0.052 * g / cm3, ncomponents = 4);
  rohacell_foam_51->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 8);
  rohacell_foam_51->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 11);
  rohacell_foam_51->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 2);
  rohacell_foam_51->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), 1);

  // making Carbon PEEK : 30 - 70 Vf.
  // https://www.quantum-polymers.com/wp-content/uploads/2017/03/QuantaPEEK-CF30.pdf
  G4Material *peek = new G4Material("PEEK", density = 1.32 * g / cm3, ncomponents = 3);
  peek->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 19);
  peek->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 12);
  peek->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 3);

  G4Material *cf = new G4Material("CF", density = 1.62 * g / cm3, ncomponents = 1);
  cf->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 1.);

  G4Material *cf30_peek70 = new G4Material("CF30_PEEK70", density = (1.32 * 0.7 + 1.62 * 0.3) * g / cm3, ncomponents = 2);
  cf30_peek70->AddMaterial(cf, fractionmass = 0.34468085);
  cf30_peek70->AddMaterial(peek, fractionmass = 0.65531915);

  // that seems to be the composition of 304 Stainless steel
  G4Material *StainlessSteel =
      new G4Material("SS304", density = 7.9 * g / cm3, ncomponents = 8);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.70105);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Cr"), 0.18);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ni"), 0.09);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.02);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.0007);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("S"), 0.0003);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.0075);
  StainlessSteel->AddElement(G4NistManager::Instance()->FindOrBuildElement("P"), 0.00045);

  G4Material *SS310 =
      new G4Material("SS310", density = 8.0 * g / cm3, ncomponents = 8);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.50455);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("Cr"), 0.25);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ni"), 0.20);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.02);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.0025);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("S"), 0.015);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.0075);
  SS310->AddElement(G4NistManager::Instance()->FindOrBuildElement("P"), 0.00045);

  // SS316 from https://www.azom.com
  G4Material *SS316 =
      new G4Material("SS316", density = 8.0 * g / cm3, ncomponents = 9);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.68095);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Cr"), 0.16);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ni"), 0.11);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.02);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mo"), 0.02);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.0008);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("S"), 0.0003);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.0075);
  SS316->AddElement(G4NistManager::Instance()->FindOrBuildElement("P"), 0.00045);

  G4Material *Steel =
      new G4Material("Steel", density = 7.86 * g / cm3, ncomponents = 5);
  Steel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.9834);
  Steel->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.014);
  Steel->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.0017);
  Steel->AddElement(G4NistManager::Instance()->FindOrBuildElement("S"), 0.00045);
  Steel->AddElement(G4NistManager::Instance()->FindOrBuildElement("P"), 0.00045);

  // a36 steel from http://www.matweb.com
  G4Material *a36 = new G4Material("Steel_A36", density = 7.85 * g / cm3, ncomponents = 5);
  a36->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.9824);
  a36->AddElement(G4NistManager::Instance()->FindOrBuildElement("Cu"), 0.002);
  a36->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.0025);
  a36->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.0103);
  a36->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.0028);

  // 1006 steel from http://www.matweb.com
  G4Material *steel_1006 = new G4Material("Steel_1006", density = 7.872 * g / cm3, ncomponents = 2);
  steel_1006->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.996);
  steel_1006->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.004);

  // from www.aalco.co.uk
  G4Material *Al5083 = new G4Material("Al5083", density = 2.65 * g / cm3, ncomponents = 3);
  Al5083->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mn"), 0.004);
  Al5083->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mg"), 0.04);
  Al5083->AddElement(G4NistManager::Instance()->FindOrBuildElement("Al"), 0.956);

  // Al 4046 from http://www.matweb.com
  G4Material *Al4046 = new G4Material("Al4046", density = 2.66 * g / cm3, ncomponents = 3);
  Al4046->AddElement(G4NistManager::Instance()->FindOrBuildElement("Al"), 0.897);
  Al4046->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.1);
  Al4046->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mg"), 0.003);

  // Al 6061T6 from http://www.matweb.com
  G4Material *Al6061T6 = new G4Material("Al6061T6", density = 2.70 * g / cm3, ncomponents = 4);
  Al6061T6->AddElement(G4NistManager::Instance()->FindOrBuildElement("Al"), 0.975);
  Al6061T6->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), 0.008);
  Al6061T6->AddElement(G4NistManager::Instance()->FindOrBuildElement("Mg"), 0.01);
  Al6061T6->AddElement(G4NistManager::Instance()->FindOrBuildElement("Fe"), 0.007);

  // E864 Pb-Scifi calorimeter
  // E864 Calorimeter is 99% Pb, 1% Antimony
  // Nuclear Instruments and Methods in Physics Research A 406 (1998) 227 258
  G4double density_e864 = (0.99 * 11.34 + 0.01 * 6.697) * g / cm3;
  G4Material *absorber_e864 = new G4Material("E864_Absorber", density_e864, 2);
  absorber_e864->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), 0.99);
  absorber_e864->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Sb"), 0.01);

  G4Material *FPC = new G4Material("FPC", 1.542 * g / cm3, 2);
  FPC->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu"), 0.0162);
  FPC->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON"), 0.9838);

  // This is an approximation for the W saturated epoxy of the EMCal.
  G4Material *W_Epoxy = new G4Material("W_Epoxy", density = 10.2 * g / cm3, ncomponents = 2);
  W_Epoxy->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_W"), fractionmass = 0.5);
  W_Epoxy->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE"), fractionmass = 0.5);

  // from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  // Epoxy (for FR4 )
  // density = 1.2*g/cm3;
  G4Material *Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, ncomponents = 2);
  Epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), natoms = 2);
  Epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 2);

  // FR4 (Glass + Epoxy)
  density = 1.86 * g / cm3;
  G4Material *FR4 = new G4Material("FR4", density, ncomponents = 2);
  FR4->AddMaterial(quartz, fractionmass = 0.528);
  FR4->AddMaterial(Epoxy, fractionmass = 0.472);
  // NOMEX (HoneyComb)
  // density from http://www.fibreglast.com/product/Nomex_Honeycomb_1562/Vacuum_Bagging_Sandwich_Core
  // 1562: 29 kg/m^3 <-- I guess it is this one
  // 2562: 48 kg/m^3
  // chemical composition http://ww2.unime.it/cdlchimind/adm/inviofile/uploads/HP_Pols2.b.pdf
  density = 29 * kg / m3;
  G4Material *NOMEX = new G4Material("NOMEX", density, ncomponents = 4);
  NOMEX->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 14);
  NOMEX->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), natoms = 10);
  NOMEX->AddElement(G4NistManager::Instance()->FindOrBuildElement("N"), natoms = 2);
  NOMEX->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), natoms = 2);
  // spacal material. Source : EICROOT/A. Kiselev
  /*
  WEpoxyMix          3  12.011 1.008 183.85  6.  1.  74.  12.18  0.029 0.002 0.969
         1  1  30.  .00001
                     0
                     */
  G4Material *Spacal_W_Epoxy =
      new G4Material("Spacal_W_Epoxy", density = 12.18 * g / cm3, ncomponents = 3);
  Spacal_W_Epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 0.029);
  Spacal_W_Epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 0.002);
  Spacal_W_Epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("W"), 0.969);
  /*
PMMA      -3  12.01 1.008 15.99  6.  1.  8.  1.19  3.6  5.7  1.4
       1  1  20.  .00001
                   0
                     */
  G4Material *PMMA =
      new G4Material("PMMA", density = 1.19 * g / cm3, ncomponents = 3);
  PMMA->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 3.6 / (3.6 + 5.7 + 1.4));
  PMMA->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 5.7 / (3.6 + 5.7 + 1.4));
  PMMA->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 1.4 / (3.6 + 5.7 + 1.4));

  // scintillator for HCal, use a new name in order to change the Birks' constant
  G4Material *Uniplast_scintillator = new G4Material("Uniplast_scintillator", 1.06 * g / cm3, ncomponents = 1);
  Uniplast_scintillator->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE"), fractionmass = 1.);

  G4Material *G10 =
      new G4Material("G10", density = 1.700 * g / cm3, ncomponents = 4);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), natoms = 1);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), natoms = 2);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 3);
  G10->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), natoms = 3);

  G4Material *CsI =
      new G4Material("CsI", density = 4.534 * g / cm3, ncomponents = 2);
  CsI->AddElement(G4NistManager::Instance()->FindOrBuildElement("Cs"), natoms = 1);
  CsI->AddElement(G4NistManager::Instance()->FindOrBuildElement("I"), natoms = 1);
  CsI->GetIonisation()->SetMeanExcitationEnergy(553.1 * eV);

  G4Material *C4F10 =
      new G4Material("C4F10", density = 0.00973 * g / cm3, ncomponents = 2);
  C4F10->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 4);
  C4F10->AddElement(G4NistManager::Instance()->FindOrBuildElement("F"), natoms = 10);

  G4Material *CF4 = new G4Material("CF4", density = 3.72 * mg / cm3, ncomponents = 2, kStateGas, 288.15 * kelvin, 1 * atmosphere);
  CF4->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), natoms = 1);
  CF4->AddElement(G4NistManager::Instance()->FindOrBuildElement("F"), natoms = 4);

  G4Element *elLu = new G4Element(name = "Lutetium", symbol = "Lu", z = 71., a = 174.97 * g / mole);
  G4Material *LSO = new G4Material("LSO",                    // its name
                                   density = 7.4 * g / cm3,  // its density
                                   ncomponents = 3);         // number of components

  LSO->AddElement(G4NistManager::Instance()->FindOrBuildElement("Si"), natoms = 1);
  LSO->AddElement(elLu, natoms = 2);
  LSO->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), natoms = 5);

  // Silver epoxy glue LOCTITE ABLESTIK 2902 for the silicon sensors and FPHX chips of INTT
  G4Material *SilverEpoxyGlue_INTT = new G4Material("SilverEpoxyGlue_INTT", density = 3.2 * g / cm3, ncomponents = 2);
  SilverEpoxyGlue_INTT->AddMaterial(Epoxy, fractionmass = 0.79);
  SilverEpoxyGlue_INTT->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Ag"), fractionmass = 0.21);

  // this here is very close but makes more sense since it uses Ne and CF4
  double G4_Ne_frac = 0.5;
  double CF4_frac = 0.5;
  const double den_G4_Ne = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne")->GetDensity();
  const double den_CF4_2 = CF4->GetDensity();
  const double den_sphenix_tpc_gas = den_G4_Ne * G4_Ne_frac + den_CF4_2 * CF4_frac;
  G4Material *sPHENIX_tpc_gas = new G4Material("sPHENIX_TPC_Gas", den_sphenix_tpc_gas, ncomponents = 2, kStateGas);
  sPHENIX_tpc_gas->AddMaterial(CF4, den_CF4_2 * CF4_frac / den_sphenix_tpc_gas);
  sPHENIX_tpc_gas->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne"), den_G4_Ne * G4_Ne_frac / den_sphenix_tpc_gas);

  // Due to supply issues, we are now expecting to use Ar CF4.
  // The fractions are tuned to produce very similar drift speed
  // and other parameters as the original NeCF4 mixture.
  double alt_G4_Ar_frac = 0.6;
  double alt_CF4_frac = 0.4;
  const double alt_den_G4_Ar = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar")->GetDensity();
  const double alt_den_CF4 = CF4->GetDensity();
  const double alt_den_sphenix_tpc_gas = alt_den_G4_Ar * alt_G4_Ar_frac + alt_den_CF4 * alt_CF4_frac;
  G4Material *alt_sPHENIX_tpc_gas = new G4Material("sPHENIX_TPC_Gas_ArCF4", alt_den_sphenix_tpc_gas, ncomponents = 2, kStateGas);
  alt_sPHENIX_tpc_gas->AddMaterial(CF4, alt_den_CF4 * alt_CF4_frac / alt_den_sphenix_tpc_gas);
  alt_sPHENIX_tpc_gas->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar"), alt_den_G4_Ar * alt_G4_Ar_frac / alt_den_sphenix_tpc_gas);

  // define P10 Gas which will be used for TPC Benchmarking
  G4Material *P10 =
      new G4Material("P10", density = 1.74 * mg / cm3, ncomponents = 3);  // @ 0K, 1atm
  P10->AddElement(G4NistManager::Instance()->FindOrBuildElement("Ar"), fractionmass = 0.9222);
  P10->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), fractionmass = 0.0623);
  P10->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), fractionmass = 0.0155);
}

void PHG4Reco::DefineRegions()
{
  // the PAI model does not work anymore in G4 10.06
  //   const G4RegionStore *theRegionStore = G4RegionStore::GetInstance();
  //   G4ProductionCuts *gcuts = new G4ProductionCuts(*(theRegionStore->GetRegion("DefaultRegionForTheWorld")->GetProductionCuts()));
  //   G4Region *tpcregion = new G4Region("REGION_TPCGAS");
  //   tpcregion->SetProductionCuts(gcuts);
  // #if G4VERSION_NUMBER >= 1033
  //   // Use this from the new G4 version 10.03 on
  //   // was commented out, crashes in 10.06 I think
  //   // add the PAI model to the TPCGAS region
  //   // undocumented, painfully digged out with debugger by tracing what
  //   // is done for command "/process/em/AddPAIRegion all TPCGAS PAI"
  // //  G4EmParameters *g4emparams = G4EmParameters::Instance();
  // //  g4emparams->AddPAIModel("all", "REGION_TPCGAS", "PAI");
  // #endif
  return;
}

PHG4Subsystem *
PHG4Reco::getSubsystem(const std::string &name)
{
  for (PHG4Subsystem *subsys : m_SubsystemList)
  {
    if (subsys->Name() == name)
    {
      if (Verbosity() > 0)
      {
        std::cout << "Found Subsystem " << name << std::endl;
      }
      return subsys;
    }
  }
  std::cout << "Could not find Subsystem " << name << std::endl;
  return nullptr;
}

void PHG4Reco::G4Verbosity(const int i)
{
  if (m_RunManager)
  {
    m_RunManager->SetVerboseLevel(i);
  }
}

void PHG4Reco::ApplyDisplayAction()
{
  if (!m_Detector)
  {
    return;
  }
  G4VPhysicalVolume *physworld = m_Detector->GetPhysicalVolume();
  m_DisplayAction->ApplyDisplayAction(physworld);
  for (PHG4Subsystem *g4sub : m_SubsystemList)
  {
    PHG4DisplayAction *action = g4sub->GetDisplayAction();
    if (action)
    {
      if (Verbosity() > 1)
      {
        std::cout << "Adding steppingaction for " << g4sub->Name() << std::endl;
      }
      action->ApplyDisplayAction(physworld);
    }
  }
}
