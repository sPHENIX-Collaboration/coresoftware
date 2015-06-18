#include "PHG4Reco.h"

#include "PHG4PrimaryGeneratorAction.h"
#include "G4TBMagneticFieldSetup.hh"
#include "PHG4PhenixDetector.h"
#include "PHG4PhenixSteppingAction.h"
#include "PHG4PhenixTrackingAction.h"
#include "PHG4PhenixEventAction.h"
#include "PHG4Subsystem.h"
#include "PHG4InEvent.h"
#include "PHG4Utils.h"

#include <g4decayer/P6DExtDecayerPhysics.hh>
#include <g4decayer/EDecayType.hh>

#include <fun4all/getClass.h>
#include <fun4all/recoConsts.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <TThread.h>

#include <CLHEP/Random/Random.h>

#include <Geant4/G4RunManager.hh>

#include <Geant4/G4VisExecutive.hh>
#include <Geant4/G4OpenGLImmediateX.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4UIExecutive.hh>
#include <Geant4/G4UImanager.hh>

#include <Geant4/G4Cerenkov.hh>
#include <Geant4/G4Scintillation.hh>
#include <Geant4/G4OpAbsorption.hh>
#include <Geant4/G4OpRayleigh.hh>
#include <Geant4/G4OpMieHG.hh>
#include <Geant4/G4OpBoundaryProcess.hh>
#include <Geant4/G4LossTableManager.hh>
#include <Geant4/G4EmSaturation.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ParticleTypes.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4OpticalPhoton.hh>
#include <Geant4/G4OpticalPhysics.hh>
#include <Geant4/G4OpWLS.hh>
#include <Geant4/G4PEEffectFluoModel.hh>
#include <Geant4/G4EmProcessOptions.hh>
#include <Geant4/G4HadronicProcessStore.hh>
#include <Geant4/G4LossTableManager.hh>

#include <Geant4/globals.hh>

#include <Geant4/G4Version.hh>

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

#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>

#include <cstdlib>
#include <cstring>
#include <cassert>

using namespace std;

static TThread *gui_thread = NULL;

// for the G4 cmd line interface
G4UImanager* UImanager = NULL;

// the gui thread
void g4guithread(void *ptr);

//_________________________________________________________________
PHG4Reco::PHG4Reco( const string &name ) :
  SubsysReco( name ),
  magfield(2),
  field_(0),
  runManager_(0),
  detector_(0),
  steppingAction_(0),
  trackingAction_(NULL),
  generatorAction_(NULL),
  visManager(0),
  _eta_coverage(1.0),
  fieldmapfile("NONE"),
  worldshape("G4Tubs"),
  worldmaterial("G4_AIR"),
  physicslist("QGSP_BERT"),
  active_decayer_(true),
  active_force_decay_(false),
  force_decay_type_(kAll),
  _timer( PHTimeServer::get()->insert_new( name ) )
{
  for (int i = 0; i < 3; i++)
    {
      WorldSize[i] = 1000.;
    }
  return;
}

//_________________________________________________________________
PHG4Reco::~PHG4Reco( void )
{
  // one can delete null pointer (it results in a nop), so checking if
  // they are non zero is not needed
  delete gui_thread;
  delete field_;
  delete runManager_;
  delete visManager;
}

//_________________________________________________________________
int PHG4Reco::Init( PHCompositeNode* topNode )
{
  if (verbosity >= 0) {
    cout << "========================= PHG4Reco::Init() ================================" << endl;
  }
  
  if (verbosity > 0) cout << "PHG4Reco::Init" << endl;
  // create GEANT run manager
  if (verbosity > 0) cout << "PHG4Reco::Init - create run manager" << endl;
  runManager_ = new G4RunManager();

  DefineMaterials();

  //setup the constant field
  if (verbosity > 0) cout << "PHG4Reco::Init - create magnetic field setup" << endl;
  if (fieldmapfile != "NONE")
    {
      field_ = new G4TBMagneticFieldSetup(fieldmapfile, mapdim) ;
      magfield = field_->get_magfield_at_000(2); // get the z coordinate at 0/0/0
      if (verbosity > 0)
	{
	  cout << "magfield in PHG4Reco: " << magfield << endl;
	}
    }
  else
    {
      field_ = new G4TBMagneticFieldSetup(magfield) ;
    }

  // create physics processes
  G4VModularPhysicsList *myphysicslist = NULL;
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
      setenv("AllowForHeavyElements","1",1);
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
       setenv("AllowForHeavyElements","1",1);
       myphysicslist = new FTFP_BERT_HP(verbosity);
     }
#endif
#ifdef HAVE_QGSP_BERT_HP
   else if (physicslist == "QGSP_BERT_HP")
     {
       setenv("AllowForHeavyElements","1",1);
       myphysicslist = new QGSP_BERT_HP(verbosity);
     }
#endif
  else
    {
      cout << "Physics List " << physicslist << " not implemented" << endl;
      gSystem->Exit(1);
    }

  if (active_decayer_) {
    P6DExtDecayerPhysics *decayer = new P6DExtDecayerPhysics();
    if (active_force_decay_) decayer->SetForceDecay(force_decay_type_);
    myphysicslist->RegisterPhysics(decayer);
  }
  runManager_->SetUserInitialization(myphysicslist);

  // initialize registered subsystems
  BOOST_FOREACH(SubsysReco * reco, subsystems_)
    {
      reco->Init( topNode );
    }
  recoConsts *rc = recoConsts::instance();

  // initialize registered subsystems
  BOOST_FOREACH(SubsysReco * reco, subsystems_)
    {
      reco->InitRun( topNode );
    }

  // create phenix detector, add subsystems, and register to GEANT
  if (verbosity > 0) cout << "PHG4Reco::Init - create detector" << endl;
  detector_ = new PHG4PhenixDetector();
  detector_->Verbosity(verbosity);
  detector_->SetWorldSizeX(WorldSize[0]*cm);
  detector_->SetWorldSizeY(WorldSize[1]*cm);
  detector_->SetWorldSizeZ(WorldSize[2]*cm);
  detector_->SetWorldShape(worldshape);
  detector_->SetWorldMaterial(worldmaterial);

  rc->set_FloatFlag("WorldSizex", WorldSize[0]);
  rc->set_FloatFlag("WorldSizey", WorldSize[1]);
  rc->set_FloatFlag("WorldSizez", WorldSize[2]);

  BOOST_FOREACH( PHG4Subsystem * g4sub, subsystems_)
    {
      detector_->AddDetector( g4sub->GetDetector() );
    }
  runManager_->SetUserInitialization( detector_ );

  setupInputEventNodeReader(topNode);
  // create main event action, add subsystemts and register to GEANT
  eventAction_ = new PHG4PhenixEventAction();

  BOOST_FOREACH( PHG4Subsystem * g4sub, subsystems_)
    {
      PHG4EventAction *evtact = g4sub->GetEventAction();
      if (evtact)
	{
	  eventAction_->AddAction(evtact);
	}
    }
  runManager_->SetUserAction(eventAction_ );

  // create main stepping action, add subsystems and register to GEANT
  steppingAction_ = new PHG4PhenixSteppingAction();
  BOOST_FOREACH( PHG4Subsystem * g4sub, subsystems_)
    {
      PHG4SteppingAction *action = g4sub->GetSteppingAction();
      if (action)
	{
	  if (verbosity > 0)
	    {
	      cout << "Adding steppingaction for " << g4sub->Name() << endl;
	    }
	  steppingAction_->AddAction( g4sub->GetSteppingAction() );
	}
    }
  runManager_->SetUserAction(steppingAction_ );

  // create main tracking action, add subsystems and register to GEANT
  trackingAction_ = new PHG4PhenixTrackingAction();
  BOOST_FOREACH( PHG4Subsystem * g4sub, subsystems_)
    {
      trackingAction_->AddAction( g4sub->GetTrackingAction() );
    }
  runManager_->SetUserAction(trackingAction_ );

  // initialize
  runManager_->Initialize();

  // add cerenkov and optical photon processes
  // cout << endl << "Ignore the next message - we implemented this correctly" << endl;
  G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
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

  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* _theParticleIterator;
  _theParticleIterator = theParticleTable->GetIterator();
  _theParticleIterator->reset();
  while( (*_theParticleIterator)() )
    {
      G4ParticleDefinition* particle = _theParticleIterator->value();
      G4String particleName = particle->GetParticleName();
      G4ProcessManager* pmanager = particle->GetProcessManager();
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
  G4ProcessManager* pmanager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
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
  //   if( G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager() ){
  //     trackingManager->SetStoreTrajectory( true );
  //   }

  // quiet some G4 print-outs (EM and Hadronic settings during first event)
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4LossTableManager::Instance()->SetVerbose(0);

  if (verbosity >= 0) {
    cout << "===========================================================================" << endl;
  }
  
  return 0;
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

int
PHG4Reco::InitUImanager()
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
int PHG4Reco::process_event( PHCompositeNode* topNode )
{
  TThread::Lock();
  // make sure Actions and subsystems have the relevant pointers set
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  generatorAction_->SetInEvent(ineve);

  BOOST_FOREACH(SubsysReco * reco, subsystems_)
    {
      if (Verbosity() >= 2)
	cout << "PHG4Reco::process_event - " << reco->Name() << "->process_event" << endl;

      try
	{
	  reco->process_event( topNode );

	}
      catch (const exception& e)
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
      cout << " PHG4Reco::process_event - " << "run one event :" << endl;
      ineve->identify();
    }
  runManager_->BeamOn( 1 );
  _timer.get()->stop();

  BOOST_FOREACH( PHG4Subsystem * g4sub, subsystems_)
    {
      if (Verbosity() >= 2)
	cout << " PHG4Reco::process_event - " << g4sub->Name() << "->process_after_geant" << endl;
      try
	{
	  g4sub->process_after_geant(topNode );
	}
      catch (const exception& e)
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


int
PHG4Reco::ResetEvent(PHCompositeNode *topNode)
{
  BOOST_FOREACH(SubsysReco * reco, subsystems_)
    {
      reco->ResetEvent( topNode );
    }
  return 0;
}

//_________________________________________________________________
int PHG4Reco::End( PHCompositeNode* )
{
  return 0;
}

int
PHG4Reco::setupInputEventNodeReader(PHCompositeNode *topNode)
{
  PHG4InEvent *ineve = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!ineve)
    {
      PHNodeIterator iter( topNode );
      PHCompositeNode *dstNode;
      dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

      ineve = new PHG4InEvent();
      PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(ineve, "PHG4INEVENT", "PHObject");
      dstNode->addNode(newNode);
    }
  generatorAction_ = new PHG4PrimaryGeneratorAction();
  runManager_->SetUserAction(generatorAction_ );
  return 0;
}

void
PHG4Reco::set_rapidity_coverage(const double eta)
{
  _eta_coverage = eta;
  PHG4Utils::SetPseudoRapidityCoverage(eta);
}

void
PHG4Reco::G4Seed(const int i)
{
  CLHEP::HepRandom::setTheSeed(i);
  return;
}

void g4guithread(void *ptr)
{
  TThread::Lock();
  G4UIExecutive* ui = new G4UIExecutive(0, NULL);
  if (ui->IsGUI() && boost::filesystem::exists( "gui.mac" ))
    {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
  TThread::UnLock();
  ui->SessionStart();
  TThread::Lock();
  delete ui;
  TThread::UnLock();
  gui_thread = NULL;
  return;
}

//____________________________________________________________________________
void
PHG4Reco::DefineMaterials()
{
  G4String symbol;             //a=mass of a mole;
  G4double density;      //z=mean number of protons;
  G4double fractionmass;

  G4int ncomponents, natoms;
  // this is for FTFP_BERT_HP where the neutron code barfs
  // if the z difference to the last known element (U) is too large
  set<G4String> ignoremat;
  ignoremat.insert("G4_Cf");

  //load all Materials from the nist database
  G4NistManager * nist = G4NistManager::Instance();
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
  //Elements needed
  G4Element *Ar  = nist->FindOrBuildElement("Ar");
  G4Element *C  = nist->FindOrBuildElement("C");
  G4Element *Cr  = nist->FindOrBuildElement("Cr");
  G4Element *Cs  = nist->FindOrBuildElement("Cs");
  G4Element *F  = nist->FindOrBuildElement("F");
  G4Element *Fe  = nist->FindOrBuildElement("Fe");
  G4Element *H  = nist->FindOrBuildElement("H");
  G4Element *I  = nist->FindOrBuildElement("I");
  G4Element *Li  = nist->FindOrBuildElement("Li");
  G4Element *Mn  = nist->FindOrBuildElement("Mn");
  G4Element *Ni  = nist->FindOrBuildElement("Ni");
  G4Element *O  = nist->FindOrBuildElement("O");
  G4Element *P  = nist->FindOrBuildElement("P");
  G4Element *S  = nist->FindOrBuildElement("S");
  G4Element *Si  = nist->FindOrBuildElement("Si");
 // making quartz
  G4Material* quartz = new G4Material
    ("Quartz", density=2.200*g/cm3, ncomponents=2);
  quartz->AddElement(Si, 1);
  quartz->AddElement(O , 2);

  // gas mixture for the MuID in fsPHENIX. CLS 02-25-14
  G4Material* IsoButane = 
    new G4Material("Isobutane", 0.00265*g/cm3, 2);
  IsoButane->AddElement( C, 4);
  IsoButane->AddElement( H, 10);

  G4Material* CO2 = new G4Material("CO2", density = 1.977e-3 * g / cm3, ncomponents=2);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);

  G4Material* MuIDgas = new G4Material("MuIDgas", density = (1.977e-3*0.92+0.00265*0.08) * g/cm3, ncomponents=2);
  MuIDgas->AddMaterial(IsoButane, fractionmass=0.08);
  MuIDgas->AddMaterial(CO2, fractionmass=0.92);
  //----

  G4Material* Sci =
    new G4Material("Scintillator", density = 1.032 * g / cm3, ncomponents = 2);
  Sci->AddElement(C, natoms = 9);
  Sci->AddElement(H, natoms = 10);

  // that seems to be the composition of 304 Stainless steel
  G4Material * StainlessSteel =
    new G4Material("SS304", density = 7.9 * g / cm3, ncomponents = 8);
  StainlessSteel->AddElement(Fe, 0.70105);
  StainlessSteel->AddElement(Cr, 0.18);
  StainlessSteel->AddElement(Ni, 0.09);
  StainlessSteel->AddElement(Mn, 0.02);
  StainlessSteel->AddElement(C, 0.0007);
  StainlessSteel->AddElement(S, 0.0003);
  StainlessSteel->AddElement(Si, 0.0075);
  StainlessSteel->AddElement(P, 0.00045);

  G4Material * SS310 =
    new G4Material("SS310", density = 8.0 * g / cm3, ncomponents = 8);
  StainlessSteel->AddElement(Fe, 0.50455);
  StainlessSteel->AddElement(Cr, 0.25);
  StainlessSteel->AddElement(Ni, 0.20);
  StainlessSteel->AddElement(Mn, 0.02);
  StainlessSteel->AddElement(C, 0.0025);
  StainlessSteel->AddElement(S, 0.015);
  StainlessSteel->AddElement(Si, 0.0075);
  StainlessSteel->AddElement(P, 0.00045);

  G4Material * Steel =
    new G4Material("Steel", density = 7.86 * g / cm3, ncomponents = 5);
  Steel->AddElement(Fe, 0.9834);
  Steel->AddElement(Mn, 0.014);
  Steel->AddElement(C, 0.0017);
  Steel->AddElement(S, 0.00045);
  Steel->AddElement(P, 0.00045);

  // This is an approximation for the W saturated epoxy of the EMCal.
  G4Material *W = nist->FindOrBuildMaterial("G4_W");
  G4Material *Epoxy = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material *W_Epoxy =
    new G4Material("W_Epoxy", density = 10.2 * g / cm3, ncomponents = 2);
  W_Epoxy->AddMaterial(W, fractionmass = 0.5);
  W_Epoxy->AddMaterial(Epoxy, fractionmass = 0.5);

  // spacal material. Source : EICROOT/A. Kiselev
  /*
  WEpoxyMix          3  12.011 1.008 183.85  6.  1.  74.  12.18  0.029 0.002 0.969
         1  1  30.  .00001
                     0
                     */
  G4Material *Spacal_W_Epoxy =
    new G4Material("Spacal_W_Epoxy", density = 12.18 * g / cm3, ncomponents = 3);
  Spacal_W_Epoxy->AddElement(nist->FindOrBuildElement("C"), 0.029);
  Spacal_W_Epoxy->AddElement(nist->FindOrBuildElement("H"), 0.002);
  Spacal_W_Epoxy->AddElement(nist->FindOrBuildElement("W"), 0.969);
  /*
PMMA      -3  12.01 1.008 15.99  6.  1.  8.  1.19  3.6  5.7  1.4
       1  1  20.  .00001
                   0
                     */
  G4Material *PMMA =
    new G4Material("PMMA", density = 1.19 * g / cm3, ncomponents = 3);
  PMMA->AddElement(nist->FindOrBuildElement("C"), 3.6/(3.6+5.7+1.4));
  PMMA->AddElement(nist->FindOrBuildElement("H"), 5.7/(3.6+5.7+1.4));
  PMMA->AddElement(nist->FindOrBuildElement("O"), 1.4/(3.6+5.7+1.4));


  G4Material* G10 =
    new G4Material("G10", density = 1.700 * g / cm3, ncomponents = 4);
  G10->AddElement(Si, natoms = 1);
  G10->AddElement(O , natoms = 2);
  G10->AddElement(C , natoms = 3);
  G10->AddElement(H , natoms = 3);

  G4Material* CsI =
    new G4Material("CsI", density = 4.534 * g / cm3, ncomponents = 2);
  CsI->AddElement(Cs, natoms = 1);
  CsI->AddElement(I , natoms = 1);
  CsI->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);

  G4Material *C4F10 =
    new G4Material("C4F10", density = 0.00973 * g / cm3, ncomponents = 2);
  C4F10->AddElement(C , natoms = 4);
  C4F10->AddElement(F , natoms = 10);

  G4Material *CF4 = new G4Material("CF4", density = 3.72 * mg / cm3, ncomponents = 2, kStateGas, 288.15 * kelvin, 1 * atmosphere);
  CF4->AddElement(C, natoms = 1);
  CF4->AddElement(F, natoms = 4);

  // radiation length Al = 8.897cm, magnet diameter 20cm
  // this is 1 radlen of al spread over the magnet by
  // weighting its density
  G4double z;
  G4double a;
  G4Material *almag =
    new G4Material("AL_MAG", z = 13., a = 26.98 * g / mole, density = (2.7 * g / cm3) * 8.897 / 20);
  // Babar magnet thickness, still 1X0
  G4Material *albabarmag =
    new G4Material("AL_BABAR_MAG", z = 13., a = 26.98 * g / mole, density = (2.7 * g / cm3) * 8.897 / 33);
  if (verbosity > 1)
    {
      cout << "adding " << almag->GetName() << endl;
      cout << "adding " << albabarmag->GetName() << endl;
    }

    {
      //! ePHENIX TPC - Jin Huang <jhuang@bnl.gov>
      //! Ref: B. Yu et al. A gem based tpc for the legs experiment. In Nuclear Science Symposium
      //! Conference Record, 2005 IEEE, volume 2, pages 924-928, 2005. doi:10.1109/NSSMIC.2005.1596405.
      G4Material * G4_Ar = nist->FindOrBuildMaterial("G4_Ar");
      G4Material * G4_CARBON_DIOXIDE = nist->FindOrBuildMaterial(
      "G4_CARBON_DIOXIDE");

      assert(G4_Ar);
      assert(G4_CARBON_DIOXIDE);
      assert(CF4);
      
      const double den_CF4 = CF4->GetDensity() * .1;
      const double den_G4_Ar = G4_Ar->GetDensity() * .8;
      const double den_G4_CARBON_DIOXIDE = G4_CARBON_DIOXIDE->GetDensity() * .1;
      const double den = den_CF4 + den_G4_Ar + den_G4_CARBON_DIOXIDE;
      
      G4Material *ePHEINX_TPC_Gas = new G4Material("ePHEINX_TPC_Gas", den,
						   ncomponents = 3, kStateGas);
      ePHEINX_TPC_Gas->AddMaterial(CF4, den_CF4 / den);
      ePHEINX_TPC_Gas->AddMaterial(G4_Ar, den_G4_Ar / den);
      ePHEINX_TPC_Gas->AddMaterial(G4_CARBON_DIOXIDE,
				   den_G4_CARBON_DIOXIDE / den);

      //
      // CF4
      //
      const G4int nEntries_CF4 = 50;

      G4double PhotonEnergy_CF4[nEntries_CF4] =
	{  5.5*eV,  5.6*eV,  5.7*eV,  5.8*eV,  5.9*eV,
	   6.0*eV,  6.1*eV,  6.2*eV,  6.3*eV,  6.4*eV,
	   6.5*eV,  6.6*eV,  6.7*eV,  6.8*eV,  6.9*eV,
	   7.0*eV,  7.1*eV,  7.2*eV,  7.3*eV,  7.4*eV,
	   7.5*eV,  7.6*eV,  7.7*eV,  7.8*eV,  7.9*eV,
	   8.0*eV,  8.1*eV,  8.2*eV,  8.4*eV,  8.6*eV,
	   8.8*eV,  9.0*eV,  9.2*eV,  9.4*eV,  9.6*eV,
	   9.8*eV, 10.0*eV, 10.2*eV, 10.4*eV, 10.6*eV,
	   10.8*eV, 11.0*eV, 11.2*eV, 11.3*eV, 11.4*eV,
	   11.5*eV, 11.6*eV, 11.7*eV, 11.8*eV, 11.9*eV };

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
         1.000662, 1.000667, 1.000672, 1.000677, 1.000682 };
            
      G4double Absorption_CF4[nEntries_CF4] =
	{ 1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m, 
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m, 
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m };
      
      G4MaterialPropertiesTable* MPT_CF4 = new G4MaterialPropertiesTable();
      
      MPT_CF4->AddProperty("RINDEX",PhotonEnergy_CF4,RefractiveIndex_CF4,nEntries_CF4)
        ->SetSpline(true);
      MPT_CF4->AddProperty("ABSLENGTH",PhotonEnergy_CF4, Absorption_CF4,nEntries_CF4)
        ->SetSpline(true);
      
      CF4->SetMaterialPropertiesTable(MPT_CF4);

      //
      // LiF
      //
      G4Material *g4_lif = nist->FindOrBuildMaterial("G4_LITHIUM_FLUORIDE");
      G4Material *LiF = new G4Material("LiF", density = 2.635 * g / cm3, ncomponents = 2);
      LiF->AddElement(Li, natoms = 1);
      LiF->AddElement(F, natoms = 1);

      if (verbosity > 0)
	{
	  cout << "g4_lif: " << g4_lif << endl;
	  cout << "LiF: " << LiF << endl;
	}

      const G4int nEntries_LiF = 50;

      G4double PhotonEnergy_LiF[nEntries_LiF] =
	{  5.5*eV,  5.6*eV,  5.7*eV,  5.8*eV,  5.9*eV,
	   6.0*eV,  6.1*eV,  6.2*eV,  6.3*eV,  6.4*eV,
	   6.5*eV,  6.6*eV,  6.7*eV,  6.8*eV,  6.9*eV,
	   7.0*eV,  7.1*eV,  7.2*eV,  7.3*eV,  7.4*eV,
	   7.5*eV,  7.6*eV,  7.7*eV,  7.8*eV,  7.9*eV,
	   8.0*eV,  8.1*eV,  8.2*eV,  8.4*eV,  8.6*eV,
	   8.8*eV,  9.0*eV,  9.2*eV,  9.4*eV,  9.6*eV,
	   9.8*eV, 10.0*eV, 10.2*eV, 10.4*eV, 10.6*eV,
	   10.8*eV, 11.0*eV, 11.2*eV, 11.3*eV, 11.4*eV,
	   11.5*eV, 11.6*eV, 11.7*eV, 11.8*eV, 11.9*eV };

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
         1.65363, 1.66236, 1.67159, 1.68139, 1.69178 };
            
      G4double Absorption_LiF[nEntries_LiF] =
	{ 1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m, 
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m, 
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,
          1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m,   1.e4*m };

      G4MaterialPropertiesTable* MPT_LiF = new G4MaterialPropertiesTable();
      
      MPT_LiF->AddProperty("RINDEX",PhotonEnergy_LiF, RefractiveIndex_LiF, nEntries_LiF)
        ->SetSpline(true);
      MPT_LiF->AddProperty("ABSLENGTH",PhotonEnergy_LiF, Absorption_LiF, nEntries_LiF)
        ->SetSpline(true);

      LiF->SetMaterialPropertiesTable(MPT_LiF);

      //
      // CsI
      //
      const G4int nEntries_CsI = 50;

      G4double PhotonEnergy_CsI[nEntries_CsI] =
	{  5.5*eV,  5.6*eV,  5.7*eV,  5.8*eV,  5.9*eV,
	   6.0*eV,  6.1*eV,  6.2*eV,  6.3*eV,  6.4*eV,
	   6.5*eV,  6.6*eV,  6.7*eV,  6.8*eV,  6.9*eV,
	   7.0*eV,  7.1*eV,  7.2*eV,  7.3*eV,  7.4*eV,
	   7.5*eV,  7.6*eV,  7.7*eV,  7.8*eV,  7.9*eV,
	   8.0*eV,  8.1*eV,  8.2*eV,  8.4*eV,  8.6*eV,
	   8.8*eV,  9.0*eV,  9.2*eV,  9.4*eV,  9.6*eV,
	   9.8*eV, 10.0*eV, 10.2*eV, 10.4*eV, 10.6*eV,
	   10.8*eV, 11.0*eV, 11.2*eV, 11.3*eV, 11.4*eV,
	   11.5*eV, 11.6*eV, 11.7*eV, 11.8*eV, 11.9*eV };

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
	 1., 1., 1., 1., 1. };
            
      G4double Absorption_CsI[nEntries_CsI] =
	{ 0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m, 
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m, 
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,
          0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m,   0.0000001*m };

      G4MaterialPropertiesTable *MPT_CsI = new G4MaterialPropertiesTable();

      MPT_CsI->AddProperty("RINDEX",PhotonEnergy_CsI, RefractiveIndex_CsI, nEntries_CsI)
        ->SetSpline(true);
      MPT_CsI->AddProperty("ABSLENGTH",PhotonEnergy_CsI, Absorption_CsI, nEntries_CsI)
        ->SetSpline(true);

      CsI->SetMaterialPropertiesTable(MPT_CsI);

    }

    // define P10 Gas which will be used for TPC Benchmarking
    G4Material* P10 = 
    new G4Material("P10", density= 1.74*mg/cm3, ncomponents=3); // @ 0K, 1atm
    P10->AddElement(Ar, fractionmass=0.9222);
    P10->AddElement(C,  fractionmass=0.0623);
    P10->AddElement(H,  fractionmass=0.0155);
}

