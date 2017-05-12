#include "PHG4Reco.h"

#include "G4TBMagneticFieldSetup.hh"
#include "PHG4InEvent.h"
#include "PHG4PhenixDetector.h"
#include "PHG4PhenixEventAction.h"
#include "PHG4PhenixSteppingAction.h"
#include "PHG4PhenixTrackingAction.h"
#include "PHG4PrimaryGeneratorAction.h"
#include "PHG4Subsystem.h"
#include "PHG4TrackingAction.h"
#include "PHG4UIsession.h"
#include "PHG4Utils.h"

#include <g4decayer/EDecayType.hh>
#include <g4decayer/P6DExtDecayerPhysics.hh>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <phgeom/PHGeomUtility.h>
#include <g4gdml/PHG4GDMLUtility.hh>

#include <TThread.h>

#include <CLHEP/Random/Random.h>

#include <Geant4/G4RunManager.hh>

#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4OpenGLImmediateX.hh>
#include <Geant4/G4StepLimiterPhysics.hh>
#include <Geant4/G4UIExecutive.hh>
#include <Geant4/G4UImanager.hh>
#include <Geant4/G4VisExecutive.hh>

#include <Geant4/G4Cerenkov.hh>
#include <Geant4/G4EmProcessOptions.hh>
#include <Geant4/G4EmSaturation.hh>
#include <Geant4/G4HadronicProcessStore.hh>
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4OpAbsorption.hh>
#include <Geant4/G4OpBoundaryProcess.hh>
#include <Geant4/G4OpMieHG.hh>
#include <Geant4/G4OpRayleigh.hh>
#include <Geant4/G4OpWLS.hh>
#include <Geant4/G4OpticalPhoton.hh>
#include <Geant4/G4OpticalPhysics.hh>
#include <Geant4/G4PEEffectFluoModel.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleTypes.hh>
#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4Scintillation.hh>

#include <Geant4/G4Version.hh>
#include <Geant4/globals.hh>

// physics lists
#include <Geant4/FTFP_BERT.hh>
#include <Geant4/LBE.hh>
#include <Geant4/QGSP_BERT.hh>
#include <Geant4/QGSP_BIC.hh>
#include <Geant4/QGSP_BIC_HP.hh>

#if G4VERSION_NUMBER <= 951
#define HAVE_LHEP
#include <Geant4/LHEP.hh>
#endif

#if G4VERSION_NUMBER >= 1001
#define HAVE_FTFP_BERT_HP
#define HAVE_QGSP_BERT_HP
#include <Geant4/FTFP_BERT_HP.hh>
#include <Geant4/QGSP_BERT_HP.hh>
#endif

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

static TThread *gui_thread = nullptr;

// for the G4 cmd line interface
G4UImanager *UImanager = nullptr;

// the gui thread
void g4guithread(void *ptr);

//_________________________________________________________________
PHG4Reco::PHG4Reco(const string &name)
  : SubsysReco(name)
  , magfield(2)
  , magfield_rescale(1.0)
  , field_(nullptr)
  , runManager_(nullptr)
  , uisession_(nullptr)
  , detector_(nullptr)
  , eventAction_(nullptr)
  , steppingAction_(nullptr)
  , trackingAction_(nullptr)
  , generatorAction_(nullptr)
  , visManager(nullptr)
  , _eta_coverage(1.0)
  , mapdim(0)
  , fieldmapfile("NONE")
  , worldshape("G4Tubs")
  , worldmaterial("G4_AIR")
  , physicslist("QGSP_BERT")
  , active_decayer_(true)
  , active_force_decay_(false)
  , force_decay_type_(kAll)
  , save_DST_geometry_(true)
  , _timer(PHTimeServer::get()->insert_new(name))
{
  for (int i = 0; i < 3; i++)
  {
    WorldSize[i] = 1000.;
  }
  return;
}

//_________________________________________________________________
PHG4Reco::~PHG4Reco(void)
{
  // one can delete null pointer (it results in a nop), so checking if
  // they are non zero is not needed
  delete gui_thread;
  delete field_;
  delete runManager_;
  delete uisession_;
  delete visManager;
  while (subsystems_.begin() != subsystems_.end())
  {
    delete subsystems_.back();
    subsystems_.pop_back();
  }
}

//_________________________________________________________________
int PHG4Reco::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
  {
    cout << "========================= PHG4Reco::Init() ================================" << endl;
  }
  unsigned int iseed = PHRandomSeed();
  cout << Name() << " G4 Random Seed: " << iseed << endl;
  G4Seed(iseed);  // fixed seed handled in PHRandomSeed()

  // create GEANT run manager
  if (verbosity > 1) cout << "PHG4Reco::Init - create run manager" << endl;

  // redirect GEANT verbosity to nowhere
  if (verbosity < 1)
  {
    G4UImanager *uimanager = G4UImanager::GetUIpointer();
    uisession_ = new PHG4UIsession();
    uimanager->SetSession(uisession_);
    uimanager->SetCoutDestination(uisession_);
  }

  runManager_ = new G4RunManager();

  DefineMaterials();

  //setup the constant field
  if (verbosity > 1) cout << "PHG4Reco::Init - create magnetic field setup" << endl;
  if (fieldmapfile != "NONE")
  {
    field_ = new G4TBMagneticFieldSetup(fieldmapfile, mapdim, magfield_rescale);
    magfield = field_->get_magfield_at_000(2);  // get the z coordinate at 0/0/0
    if (verbosity > 1)
    {
      cout << "magfield in PHG4Reco: " << magfield << endl;
    }
  }
  else
  {
    field_ = new G4TBMagneticFieldSetup(magfield * magfield_rescale);
  }

  // create physics processes
  G4VModularPhysicsList *myphysicslist = nullptr;
  if (physicslist == "FTFP_BERT")
  {
    myphysicslist = new FTFP_BERT(verbosity);
  }
  else if (physicslist == "QGSP_BERT")
  {
    myphysicslist = new QGSP_BERT(verbosity);
  }
  else if (physicslist == "QGSP_BIC")
  {
    myphysicslist = new QGSP_BIC(verbosity);
  }
  else if (physicslist == "QGSP_BIC_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new QGSP_BIC_HP(verbosity);
  }
#ifdef HAVE_LHEP
  else if (physicslist == "LHEP")
  {
    myphysicslist = new LHEP(verbosity);
  }
#endif
#ifdef HAVE_FTFP_BERT_HP
  else if (physicslist == "FTFP_BERT_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new FTFP_BERT_HP(verbosity);
  }
#endif
#ifdef HAVE_QGSP_BERT_HP
  else if (physicslist == "QGSP_BERT_HP")
  {
    setenv("AllowForHeavyElements", "1", 1);
    myphysicslist = new QGSP_BERT_HP(verbosity);
  }
#endif
  else
  {
    cout << "Physics List " << physicslist << " not implemented" << endl;
    gSystem->Exit(1);
  }

  if (active_decayer_)
  {
    P6DExtDecayerPhysics *decayer = new P6DExtDecayerPhysics();
    if (active_force_decay_) decayer->SetForceDecay(force_decay_type_);
    myphysicslist->RegisterPhysics(decayer);
  }

  myphysicslist->RegisterPhysics(new G4StepLimiterPhysics());
  runManager_->SetUserInitialization(myphysicslist);

  // initialize registered subsystems
  BOOST_FOREACH (SubsysReco *reco, subsystems_)
  {
    reco->Init(topNode);
  }

  if (verbosity > 0)
  {
    cout << "===========================================================================" << endl;
  }
  return 0;
}

int PHG4Reco::InitRun(PHCompositeNode *topNode)
{
  // this is a dumb protection against executing this twice.
  // we have cases (currently detector display or material scan) where
  // we need the detector bu have not run any events (who wants to wait
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
  if (verbosity > 0)
  {
    cout << "========================= PHG4Reco::InitRun() ================================" << endl;
  }

  recoConsts *rc = recoConsts::instance();

  // initialize registered subsystems
  BOOST_FOREACH (SubsysReco *reco, subsystems_)
  {
    reco->InitRun(topNode);
  }

  // create phenix detector, add subsystems, and register to GEANT
  if (verbosity > 1) cout << "PHG4Reco::Init - create detector" << endl;
  detector_ = new PHG4PhenixDetector();
  detector_->Verbosity(verbosity);
  detector_->SetWorldSizeX(WorldSize[0] * cm);
  detector_->SetWorldSizeY(WorldSize[1] * cm);
  detector_->SetWorldSizeZ(WorldSize[2] * cm);
  detector_->SetWorldShape(worldshape);
  detector_->SetWorldMaterial(worldmaterial);

  rc->set_FloatFlag("WorldSizex", WorldSize[0]);
  rc->set_FloatFlag("WorldSizey", WorldSize[1]);
  rc->set_FloatFlag("WorldSizez", WorldSize[2]);

  BOOST_FOREACH (PHG4Subsystem *g4sub, subsystems_)
  {
    detector_->AddDetector(g4sub->GetDetector());
  }
  runManager_->SetUserInitialization(detector_);

  setupInputEventNodeReader(topNode);
  // create main event action, add subsystemts and register to GEANT
  eventAction_ = new PHG4PhenixEventAction();

  BOOST_FOREACH (PHG4Subsystem *g4sub, subsystems_)
  {
    PHG4EventAction *evtact = g4sub->GetEventAction();
    if (evtact)
    {
      eventAction_->AddAction(evtact);
    }
  }
  runManager_->SetUserAction(eventAction_);

  // create main stepping action, add subsystems and register to GEANT
  steppingAction_ = new PHG4PhenixSteppingAction();
  BOOST_FOREACH (PHG4Subsystem *g4sub, subsystems_)
  {
    PHG4SteppingAction *action = g4sub->GetSteppingAction();
    if (action)
    {
      if (verbosity > 1)
      {
        cout << "Adding steppingaction for " << g4sub->Name() << endl;
      }
      steppingAction_->AddAction(g4sub->GetSteppingAction());
    }
  }
  runManager_->SetUserAction(steppingAction_);

  // create main tracking action, add subsystems and register to GEANT
  trackingAction_ = new PHG4PhenixTrackingAction();
  BOOST_FOREACH (PHG4Subsystem *g4sub, subsystems_)
  {
    trackingAction_->AddAction(g4sub->GetTrackingAction());

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

  runManager_->SetUserAction(trackingAction_);

  // initialize
  runManager_->Initialize();

  // add cerenkov and optical photon processes
  // cout << endl << "Ignore the next message - we implemented this correctly" << endl;
  G4Cerenkov *theCerenkovProcess = new G4Cerenkov("Cerenkov");
  // cout << "End of bogus warning message" << endl << endl;
  // G4Scintillation* theScintillationProcess      = new G4Scintillation("Scintillation");

  /*
    if (verbosity > 0)
    {
    // This segfaults
    theCerenkovProcess->DumpPhysicsTable();
    }
  */
  theCerenkovProcess->SetMaxNumPhotonsPerStep(100);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);

  // theScintillationProcess->SetScintillationYieldFactor(1.);
  // theScintillationProcess->SetTrackSecondariesFirst(true);

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
    // if (theScintillationProcess->IsApplicable(*particle))
    // {
    //   pmanager->AddProcess(theScintillationProcess);
    //   pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
    //   pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    // }
  }
  G4ProcessManager *pmanager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  // G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
  pmanager->AddDiscreteProcess(new G4OpAbsorption());
  pmanager->AddDiscreteProcess(new G4OpRayleigh());
  pmanager->AddDiscreteProcess(new G4OpMieHG());
  pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());
  pmanager->AddDiscreteProcess(new G4OpWLS());
  pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
  // pmanager->DumpInfo();

  // needs large amount of memory which kills central hijing events
  // store generated trajectories
  //if( G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager() ){
  //  trackingManager->SetStoreTrajectory( true );
  //}

  // quiet some G4 print-outs (EM and Hadronic settings during first event)
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4LossTableManager::Instance()->SetVerbose(0);

  if ((verbosity < 1) && (uisession_))
  {
    uisession_->Verbosity(1);  // let messages after setup come through
  }

  // Geometry export to DST
  if (save_DST_geometry_)
  {
    const string filename =
        PHGeomUtility::
            GenerateGeometryFileName("gdml");
    cout << "PHG4Reco::InitRun - export geometry to DST via tmp file " << filename << endl;

    Dump_GDML(filename);

    PHGeomUtility::ImportGeomFile(topNode, filename);

    PHGeomUtility::RemoveGeometryFile(filename);
  }

  if (verbosity > 0)
  {
    cout << "===========================================================================" << endl;
  }

  return 0;
}

//________________________________________________________________
//Dump TGeo File
void PHG4Reco::Dump_GDML(const std::string &filename)
{
  PHG4GDMLUtility :: Dump_GDML(filename , detector_->GetPhysicalVolume());
}

//_________________________________________________________________
int PHG4Reco::ApplyCommand(const std::string &cmd)
{
  InitUImanager();
  int iret = UImanager->ApplyCommand(cmd.c_str());
  return iret;
}
//_________________________________________________________________

int PHG4Reco::StartGui()
{
  if (!gui_thread)
  {
    InitUImanager();
    gui_thread = new TThread("G4Gui", g4guithread);
    gui_thread->Run();
    return 0;
  }
  return 1;
}

int PHG4Reco::InitUImanager()
{
  if (!UImanager)
  {
    // Get the pointer to the User Interface manager
    // Initialize visualization
    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager = G4UImanager::GetUIpointer();
  }
  return 0;
}

//_________________________________________________________________
int PHG4Reco::process_event(PHCompositeNode *topNode)
{
  TThread::Lock();
  // make sure Actions and subsystems have the relevant pointers set
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  generatorAction_->SetInEvent(ineve);

  BOOST_FOREACH (SubsysReco *reco, subsystems_)
  {
    if (Verbosity() >= 2)
      cout << "PHG4Reco::process_event - " << reco->Name() << "->process_event" << endl;

    try
    {
      reco->process_event(topNode);
    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during process_event from "
           << reco->Name() << endl;
      cout << "error: " << e.what() << endl;
      TThread::UnLock();
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  _timer.get()->restart();

  // run one event
  if (Verbosity() >= 2)
  {
    cout << " PHG4Reco::process_event - "
         << "run one event :" << endl;
    ineve->identify();
  }
  runManager_->BeamOn(1);
  _timer.get()->stop();

  BOOST_FOREACH (PHG4Subsystem *g4sub, subsystems_)
  {
    if (Verbosity() >= 2)
      cout << " PHG4Reco::process_event - " << g4sub->Name() << "->process_after_geant" << endl;
    try
    {
      g4sub->process_after_geant(topNode);
    }
    catch (const exception &e)
    {
      cout << PHWHERE << " caught exception thrown during process_after_geant from "
           << g4sub->Name() << endl;
      cout << "error: " << e.what() << endl;
      TThread::UnLock();
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  TThread::UnLock();
  return 0;
}

int PHG4Reco::ResetEvent(PHCompositeNode *topNode)
{
  BOOST_FOREACH (SubsysReco *reco, subsystems_)
  {
    reco->ResetEvent(topNode);
  }
  return 0;
}

//_________________________________________________________________
int PHG4Reco::End(PHCompositeNode *)
{
  return 0;
}

void PHG4Reco::Print(const std::string &what) const
{
  BOOST_FOREACH (SubsysReco *reco, subsystems_)
  {
    if (what.empty() || what == "ALL" || (reco->Name()).find(what) != string::npos)
    {
      cout << "Printing " << reco->Name() << endl;
      reco->Print(what);
      cout << "---------------------------" << endl;
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
  generatorAction_ = new PHG4PrimaryGeneratorAction();
  runManager_->SetUserAction(generatorAction_);
  return 0;
}

void PHG4Reco::setGeneratorAction(G4VUserPrimaryGeneratorAction *action)
{
  if (runManager_)
  {
    runManager_->SetUserAction(action);
  }
  return;
}

void PHG4Reco::set_rapidity_coverage(const double eta)
{
  _eta_coverage = eta;
  PHG4Utils::SetPseudoRapidityCoverage(eta);
}

void PHG4Reco::G4Seed(const unsigned int i)
{
  CLHEP::HepRandom::setTheSeed(i);
  return;
}

void g4guithread(void *ptr)
{
  TThread::Lock();
  G4UIExecutive *ui = new G4UIExecutive(0, nullptr);
  if (ui->IsGUI() && boost::filesystem::exists("gui.mac"))
  {
    UImanager->ApplyCommand("/control/execute gui.mac");
  }
  TThread::UnLock();
  ui->SessionStart();
  TThread::Lock();
  delete ui;
  TThread::UnLock();
  gui_thread = nullptr;
  return;
}

//____________________________________________________________________________
void PHG4Reco::DefineMaterials()
{
  G4String symbol;   //a=mass of a mole;
  G4double density;  //z=mean number of protons;
  G4double fractionmass;

  G4int ncomponents, natoms;
  // this is for FTFP_BERT_HP where the neutron code barfs
  // if the z difference to the last known element (U) is too large
  set<G4String> ignoremat;
  ignoremat.insert("G4_Cf");

  //load all Materials from the nist database
  G4NistManager *nist = G4NistManager::Instance();
  vector<G4String> matnames = nist->GetNistMaterialNames();
  while (matnames.begin() != matnames.end())
  {
    G4String mat = matnames.back();
    if (ignoremat.find(mat) == ignoremat.end())
    {
      nist->FindOrBuildMaterial(mat);
    }
    matnames.pop_back();
  }
  // home made compounds
  // making quartz
  G4Material *quartz = new G4Material("Quartz", density = 2.200 * g / cm3, ncomponents = 2);
  quartz->AddElement(G4Element::GetElement("Si"), 1);
  quartz->AddElement(G4Element::GetElement("O"), 2);

  // gas mixture for the MuID in fsPHENIX. CLS 02-25-14
  G4Material *IsoButane = new G4Material("Isobutane", 0.00265 * g / cm3, 2);
  IsoButane->AddElement(G4Element::GetElement("C"), 4);
  IsoButane->AddElement(G4Element::GetElement("H"), 10);

  G4Material *MuIDgas = new G4Material("MuIDgas", density = (1.977e-3 * 0.92 + 0.00265 * 0.08) * g / cm3, ncomponents = 2);
  MuIDgas->AddMaterial(IsoButane, fractionmass = 0.08);
  MuIDgas->AddMaterial(G4Material::GetMaterial("G4_CARBON_DIOXIDE"), fractionmass = 0.92);

  // that seems to be the composition of 304 Stainless steel
  G4Material *StainlessSteel =
      new G4Material("SS304", density = 7.9 * g / cm3, ncomponents = 8);
  StainlessSteel->AddElement(G4Element::GetElement("Fe"), 0.70105);
  StainlessSteel->AddElement(G4Element::GetElement("Cr"), 0.18);
  StainlessSteel->AddElement(G4Element::GetElement("Ni"), 0.09);
  StainlessSteel->AddElement(G4Element::GetElement("Mn"), 0.02);
  StainlessSteel->AddElement(G4Element::GetElement("C"), 0.0007);
  StainlessSteel->AddElement(G4Element::GetElement("S"), 0.0003);
  StainlessSteel->AddElement(G4Element::GetElement("Si"), 0.0075);
  StainlessSteel->AddElement(G4Element::GetElement("P"), 0.00045);

  G4Material *SS310 =
      new G4Material("SS310", density = 8.0 * g / cm3, ncomponents = 8);
  SS310->AddElement(G4Element::GetElement("Fe"), 0.50455);
  SS310->AddElement(G4Element::GetElement("Cr"), 0.25);
  SS310->AddElement(G4Element::GetElement("Ni"), 0.20);
  SS310->AddElement(G4Element::GetElement("Mn"), 0.02);
  SS310->AddElement(G4Element::GetElement("C"), 0.0025);
  SS310->AddElement(G4Element::GetElement("S"), 0.015);
  SS310->AddElement(G4Element::GetElement("Si"), 0.0075);
  SS310->AddElement(G4Element::GetElement("P"), 0.00045);

  G4Material *Steel =
      new G4Material("Steel", density = 7.86 * g / cm3, ncomponents = 5);
  Steel->AddElement(G4Element::GetElement("Fe"), 0.9834);
  Steel->AddElement(G4Element::GetElement("Mn"), 0.014);
  Steel->AddElement(G4Element::GetElement("C"), 0.0017);
  Steel->AddElement(G4Element::GetElement("S"), 0.00045);
  Steel->AddElement(G4Element::GetElement("P"), 0.00045);

  // a36 steel from http://www.matweb.com
  G4Material *a36 = new G4Material("Steel_A36", density = 7.85 * g / cm3, ncomponents = 5);
  a36->AddElement(G4Element::GetElement("Fe"), 0.9824);
  a36->AddElement(G4Element::GetElement("Cu"), 0.002);
  a36->AddElement(G4Element::GetElement("C"), 0.0025);
  a36->AddElement(G4Element::GetElement("Mn"), 0.0103);
  a36->AddElement(G4Element::GetElement("Si"), 0.0028);

  // 1006 steel from http://www.matweb.com
  G4Material *steel_1006 = new G4Material("Steel_1006", density = 7.872 * g / cm3, ncomponents = 2);
  steel_1006->AddElement(G4Element::GetElement("Fe"), 0.996);
  steel_1006->AddElement(G4Element::GetElement("Mn"), 0.004);

  // from www.aalco.co.uk
  G4Material *Al5083 = new G4Material("Al5083", density = 2.65 * g / cm3, ncomponents = 3);
  Al5083->AddElement(G4Element::GetElement("Mn"), 0.004);
  Al5083->AddElement(G4Element::GetElement("Mg"), 0.04);
  Al5083->AddElement(G4Element::GetElement("Al"), 0.956);

  G4Material *FPC = new G4Material("FPC", 1.542 * g / cm3, 2);
  FPC->AddMaterial(G4Material::GetMaterial("G4_Cu"), 0.0162);
  FPC->AddMaterial(G4Material::GetMaterial("G4_KAPTON"), 0.9838);

  // This is an approximation for the W saturated epoxy of the EMCal.
  G4Material *W_Epoxy = new G4Material("W_Epoxy", density = 10.2 * g / cm3, ncomponents = 2);
  W_Epoxy->AddMaterial(G4Material::GetMaterial("G4_W"), fractionmass = 0.5);
  W_Epoxy->AddMaterial(G4Material::GetMaterial("G4_POLYSTYRENE"), fractionmass = 0.5);

  //from http://www.physi.uni-heidelberg.de/~adler/TRD/TRDunterlagen/RadiatonLength/tgc2.htm
  //Epoxy (for FR4 )
  //density = 1.2*g/cm3;
  G4Material *Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, ncomponents = 2);
  Epoxy->AddElement(G4Element::GetElement("H"), natoms = 2);
  Epoxy->AddElement(G4Element::GetElement("C"), natoms = 2);

  //FR4 (Glass + Epoxy)
  density = 1.86 * g / cm3;
  G4Material *FR4 = new G4Material("FR4", density, ncomponents = 2);
  FR4->AddMaterial(quartz, fractionmass = 0.528);
  FR4->AddMaterial(Epoxy, fractionmass = 0.472);
  // spacal material. Source : EICROOT/A. Kiselev
  /*
  WEpoxyMix          3  12.011 1.008 183.85  6.  1.  74.  12.18  0.029 0.002 0.969
         1  1  30.  .00001
                     0
                     */
  G4Material *Spacal_W_Epoxy =
      new G4Material("Spacal_W_Epoxy", density = 12.18 * g / cm3, ncomponents = 3);
  Spacal_W_Epoxy->AddElement(G4Element::GetElement("C"), 0.029);
  Spacal_W_Epoxy->AddElement(G4Element::GetElement("H"), 0.002);
  Spacal_W_Epoxy->AddElement(G4Element::GetElement("W"), 0.969);
  /*
PMMA      -3  12.01 1.008 15.99  6.  1.  8.  1.19  3.6  5.7  1.4
       1  1  20.  .00001
                   0
                     */
  G4Material *PMMA =
      new G4Material("PMMA", density = 1.19 * g / cm3, ncomponents = 3);
  PMMA->AddElement(G4Element::GetElement("C"), 3.6 / (3.6 + 5.7 + 1.4));
  PMMA->AddElement(G4Element::GetElement("H"), 5.7 / (3.6 + 5.7 + 1.4));
  PMMA->AddElement(G4Element::GetElement("O"), 1.4 / (3.6 + 5.7 + 1.4));

  G4Material *G10 =
      new G4Material("G10", density = 1.700 * g / cm3, ncomponents = 4);
  G10->AddElement(G4Element::GetElement("Si"), natoms = 1);
  G10->AddElement(G4Element::GetElement("O"), natoms = 2);
  G10->AddElement(G4Element::GetElement("C"), natoms = 3);
  G10->AddElement(G4Element::GetElement("H"), natoms = 3);

  G4Material *CsI =
      new G4Material("CsI", density = 4.534 * g / cm3, ncomponents = 2);
  CsI->AddElement(G4Element::GetElement("Cs"), natoms = 1);
  CsI->AddElement(G4Element::GetElement("I"), natoms = 1);
  CsI->GetIonisation()->SetMeanExcitationEnergy(553.1 * eV);

  G4Material *C4F10 =
      new G4Material("C4F10", density = 0.00973 * g / cm3, ncomponents = 2);
  C4F10->AddElement(G4Element::GetElement("C"), natoms = 4);
  C4F10->AddElement(G4Element::GetElement("F"), natoms = 10);

  G4Material *CF4 = new G4Material("CF4", density = 3.72 * mg / cm3, ncomponents = 2, kStateGas, 288.15 * kelvin, 1 * atmosphere);
  CF4->AddElement(G4Element::GetElement("C"), natoms = 1);
  CF4->AddElement(G4Element::GetElement("F"), natoms = 4);

  //! ePHENIX TPC - Jin Huang <jhuang@bnl.gov>
  //! Ref: B. Yu et al. A gem based tpc for the legs experiment. In Nuclear Science Symposium
  //! Conference Record, 2005 IEEE, volume 2, pages 924-928, 2005. doi:10.1109/NSSMIC.2005.1596405.

  const double den_CF4 = CF4->GetDensity() * .1;
  const double den_G4_Ar = G4Material::GetMaterial("G4_Ar")->GetDensity() * .8;
  const double den_G4_CARBON_DIOXIDE = G4Material::GetMaterial("G4_CARBON_DIOXIDE")->GetDensity() * .1;
  const double den = den_CF4 + den_G4_Ar + den_G4_CARBON_DIOXIDE;

  G4Material *ePHEINX_TPC_Gas = new G4Material("ePHEINX_TPC_Gas", den,
                                               ncomponents = 3, kStateGas);
  ePHEINX_TPC_Gas->AddMaterial(CF4, den_CF4 / den);
  ePHEINX_TPC_Gas->AddMaterial(G4Material::GetMaterial("G4_Ar"), den_G4_Ar / den);
  ePHEINX_TPC_Gas->AddMaterial(G4Material::GetMaterial("G4_CARBON_DIOXIDE"),
                               den_G4_CARBON_DIOXIDE / den);
  // cross checked with original implementation made up of Ne,C,F
  // this here is very close but makes more sense since it uses Ne and CF4
  double G4_Ne_frac = 0.9;
  double CF4_frac = 0.1;
  const double den_G4_Ne = G4Material::GetMaterial("G4_Ne")->GetDensity();
  const double den_CF4_2 = CF4->GetDensity();
  const double den_sphenix_tpc_gas = den_G4_Ne * G4_Ne_frac + den_CF4_2 * CF4_frac;
  G4Material *sPHENIX_tpc_gas = new G4Material("sPHENIX_TPC_Gas", den_sphenix_tpc_gas, ncomponents = 2, kStateGas);
  sPHENIX_tpc_gas->AddMaterial(CF4, den_CF4_2 * CF4_frac / den_sphenix_tpc_gas);
  sPHENIX_tpc_gas->AddMaterial(G4Material::GetMaterial("G4_Ne"), den_G4_Ne * G4_Ne_frac / den_sphenix_tpc_gas);
  //
  // CF4
  //
  const G4int nEntries_CF4 = 50;

  G4double PhotonEnergy_CF4[nEntries_CF4] =
      {5.5 * eV, 5.6 * eV, 5.7 * eV, 5.8 * eV, 5.9 * eV,
       6.0 * eV, 6.1 * eV, 6.2 * eV, 6.3 * eV, 6.4 * eV,
       6.5 * eV, 6.6 * eV, 6.7 * eV, 6.8 * eV, 6.9 * eV,
       7.0 * eV, 7.1 * eV, 7.2 * eV, 7.3 * eV, 7.4 * eV,
       7.5 * eV, 7.6 * eV, 7.7 * eV, 7.8 * eV, 7.9 * eV,
       8.0 * eV, 8.1 * eV, 8.2 * eV, 8.4 * eV, 8.6 * eV,
       8.8 * eV, 9.0 * eV, 9.2 * eV, 9.4 * eV, 9.6 * eV,
       9.8 * eV, 10.0 * eV, 10.2 * eV, 10.4 * eV, 10.6 * eV,
       10.8 * eV, 11.0 * eV, 11.2 * eV, 11.3 * eV, 11.4 * eV,
       11.5 * eV, 11.6 * eV, 11.7 * eV, 11.8 * eV, 11.9 * eV};

  G4double RefractiveIndex_CF4[nEntries_CF4] =
      {1.000480, 1.000482, 1.000483, 1.000485, 1.000486,
       1.000488, 1.000490, 1.000491, 1.000493, 1.000495,
       1.000497, 1.000498, 1.000500, 1.000502, 1.000504,
       1.000506, 1.000508, 1.000510, 1.000512, 1.000514,
       1.000517, 1.000519, 1.000521, 1.000524, 1.000526,
       1.000529, 1.000531, 1.000534, 1.000539, 1.000545,
       1.000550, 1.000557, 1.000563, 1.000570, 1.000577,
       1.000584, 1.000592, 1.000600, 1.000608, 1.000617,
       1.000626, 1.000636, 1.000646, 1.000652, 1.000657,
       1.000662, 1.000667, 1.000672, 1.000677, 1.000682};

  G4double Absorption_CF4[nEntries_CF4] =
      {1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m};

  G4MaterialPropertiesTable *MPT_CF4 = new G4MaterialPropertiesTable();

  MPT_CF4->AddProperty("RINDEX", PhotonEnergy_CF4, RefractiveIndex_CF4, nEntries_CF4)
      ->SetSpline(true);
  MPT_CF4->AddProperty("ABSLENGTH", PhotonEnergy_CF4, Absorption_CF4, nEntries_CF4)
      ->SetSpline(true);

  CF4->SetMaterialPropertiesTable(MPT_CF4);

  //
  // LiF
  //
  G4Material *g4_lif = nist->FindOrBuildMaterial("G4_LITHIUM_FLUORIDE");
  G4Material *LiF = new G4Material("LiF", density = 2.635 * g / cm3, ncomponents = 2);
  LiF->AddElement(G4Element::GetElement("Li"), natoms = 1);
  LiF->AddElement(G4Element::GetElement("F"), natoms = 1);

  if (verbosity > 1)
  {
    cout << "g4_lif: " << g4_lif << endl;
    cout << "LiF: " << LiF << endl;
  }

  const G4int nEntries_LiF = 50;

  G4double PhotonEnergy_LiF[nEntries_LiF] =
      {5.5 * eV, 5.6 * eV, 5.7 * eV, 5.8 * eV, 5.9 * eV,
       6.0 * eV, 6.1 * eV, 6.2 * eV, 6.3 * eV, 6.4 * eV,
       6.5 * eV, 6.6 * eV, 6.7 * eV, 6.8 * eV, 6.9 * eV,
       7.0 * eV, 7.1 * eV, 7.2 * eV, 7.3 * eV, 7.4 * eV,
       7.5 * eV, 7.6 * eV, 7.7 * eV, 7.8 * eV, 7.9 * eV,
       8.0 * eV, 8.1 * eV, 8.2 * eV, 8.4 * eV, 8.6 * eV,
       8.8 * eV, 9.0 * eV, 9.2 * eV, 9.4 * eV, 9.6 * eV,
       9.8 * eV, 10.0 * eV, 10.2 * eV, 10.4 * eV, 10.6 * eV,
       10.8 * eV, 11.0 * eV, 11.2 * eV, 11.3 * eV, 11.4 * eV,
       11.5 * eV, 11.6 * eV, 11.7 * eV, 11.8 * eV, 11.9 * eV};

  G4double RefractiveIndex_LiF[nEntries_LiF] =
      {1.42709, 1.42870, 1.42998, 1.43177, 1.43368,
       1.43520, 1.43736, 1.43907, 1.44088, 1.44279,
       1.44481, 1.44694, 1.44920, 1.45161, 1.45329,
       1.45595, 1.45781, 1.46077, 1.46285, 1.46503,
       1.46849, 1.47093, 1.47349, 1.47618, 1.47901,
       1.48198, 1.48511, 1.48841, 1.49372, 1.50152,
       1.50799, 1.51509, 1.52290, 1.53152, 1.54108,
       1.54805, 1.55954, 1.56799, 1.58202, 1.59243,
       1.60382, 1.61632, 1.63010, 1.63753, 1.64536,
       1.65363, 1.66236, 1.67159, 1.68139, 1.69178};

  G4double Absorption_LiF[nEntries_LiF] =
      {1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m,
       1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m, 1.e4 * m};

  G4MaterialPropertiesTable *MPT_LiF = new G4MaterialPropertiesTable();

  MPT_LiF->AddProperty("RINDEX", PhotonEnergy_LiF, RefractiveIndex_LiF, nEntries_LiF)
      ->SetSpline(true);
  MPT_LiF->AddProperty("ABSLENGTH", PhotonEnergy_LiF, Absorption_LiF, nEntries_LiF)
      ->SetSpline(true);

  LiF->SetMaterialPropertiesTable(MPT_LiF);

  //
  // CsI
  //
  const G4int nEntries_CsI = 50;

  G4double PhotonEnergy_CsI[nEntries_CsI] =
      {5.5 * eV, 5.6 * eV, 5.7 * eV, 5.8 * eV, 5.9 * eV,
       6.0 * eV, 6.1 * eV, 6.2 * eV, 6.3 * eV, 6.4 * eV,
       6.5 * eV, 6.6 * eV, 6.7 * eV, 6.8 * eV, 6.9 * eV,
       7.0 * eV, 7.1 * eV, 7.2 * eV, 7.3 * eV, 7.4 * eV,
       7.5 * eV, 7.6 * eV, 7.7 * eV, 7.8 * eV, 7.9 * eV,
       8.0 * eV, 8.1 * eV, 8.2 * eV, 8.4 * eV, 8.6 * eV,
       8.8 * eV, 9.0 * eV, 9.2 * eV, 9.4 * eV, 9.6 * eV,
       9.8 * eV, 10.0 * eV, 10.2 * eV, 10.4 * eV, 10.6 * eV,
       10.8 * eV, 11.0 * eV, 11.2 * eV, 11.3 * eV, 11.4 * eV,
       11.5 * eV, 11.6 * eV, 11.7 * eV, 11.8 * eV, 11.9 * eV};

  G4double RefractiveIndex_CsI[nEntries_CsI] =
      {1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.};

  G4double Absorption_CsI[nEntries_CsI] =
      {0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m,
       0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m, 0.0000001 * m};

  G4MaterialPropertiesTable *MPT_CsI = new G4MaterialPropertiesTable();

  MPT_CsI->AddProperty("RINDEX", PhotonEnergy_CsI, RefractiveIndex_CsI, nEntries_CsI)
      ->SetSpline(true);
  MPT_CsI->AddProperty("ABSLENGTH", PhotonEnergy_CsI, Absorption_CsI, nEntries_CsI)
      ->SetSpline(true);

  CsI->SetMaterialPropertiesTable(MPT_CsI);

  // define P10 Gas which will be used for TPC Benchmarking
  G4Material *P10 =
      new G4Material("P10", density = 1.74 * mg / cm3, ncomponents = 3);  // @ 0K, 1atm
  P10->AddElement(G4Element::GetElement("Ar"), fractionmass = 0.9222);
  P10->AddElement(G4Element::GetElement("C"), fractionmass = 0.0623);
  P10->AddElement(G4Element::GetElement("H"), fractionmass = 0.0155);
}

PHG4Subsystem *
PHG4Reco::getSubsystem(const string &name)
{
  BOOST_FOREACH (PHG4Subsystem *subsys, subsystems_)
  {
    if (subsys->Name() == name)
    {
      if (verbosity > 0)
      {
        cout << "Found Subsystem " << name << endl;
      }
      return subsys;
    }
  }
  cout << "Could not find Subsystem " << name << endl;
  return nullptr;
}

void
PHG4Reco::G4Verbosity(const int i)
{
  if (runManager_) 
  {
    runManager_->SetVerboseLevel(i);
  }
}
