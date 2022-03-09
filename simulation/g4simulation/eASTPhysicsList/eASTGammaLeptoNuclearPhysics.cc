////////////////////////////////////////////////////////////////////////////////
//
//  eASTGammaLeptoNuclearPhysics.cc
//  Description: Gamma-nuclear, electro-nuclear and muon-nuclear physics
//                 constructor for eASTPhysicsList
//
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////


#include "eASTGammaLeptoNuclearPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "G4PhotoNuclearProcess.hh"
#else
#include "G4HadronInelasticProcess.hh"
#include "G4PhotoNuclearCrossSection.hh"
#include "G4HadronicProcess.hh"
#endif
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4MuonNuclearProcess.hh"

#include "G4CascadeInterface.hh"
#include "G4ElectroVDNuclearModel.hh"
#include "G4MuonVDNuclearModel.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSMFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"

#include "G4SystemOfUnits.hh"

#if G4VERSION_NUMBER < 1100
eASTGammaLeptoNuclearPhysics::eASTGammaLeptoNuclearPhysics()
{}

eASTGammaLeptoNuclearPhysics::~eASTGammaLeptoNuclearPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}
#else
eASTGammaLeptoNuclearPhysics::eASTGammaLeptoNuclearPhysics()
: G4VPhysicsConstructor("eASTGammaLeptoNuclear")
{;}

eASTGammaLeptoNuclearPhysics::~eASTGammaLeptoNuclearPhysics()
{;}
#endif

void eASTGammaLeptoNuclearPhysics::ConstructProcess()
{
  // Use Bertini cascade for low energies
  G4CascadeInterface* theGammaReaction = new G4CascadeInterface;
  theGammaReaction->SetMinEnergy(0.0);
  theGammaReaction->SetMaxEnergy(3.5*GeV);

  // Use QGSP for high energies
  qgsp = new G4TheoFSGenerator("QGSP");
  stringModel = new G4QGSModel<G4GammaParticipants>;
  stringDecay =
    new G4ExcitedStringDecay(fragModel = new G4QGSMFragmentation);
  stringModel->SetFragmentationModel(stringDecay);
  preCompoundModel = new G4GeneratorPrecompoundInterface();

  qgsp->SetHighEnergyGenerator(stringModel);
  qgsp->SetTransport(preCompoundModel); 
  qgsp->SetMinEnergy(3*GeV);
  qgsp->SetMaxEnergy(100*TeV);

  // Lepto-nuclear models
  G4ElectroVDNuclearModel* evdn = new G4ElectroVDNuclearModel;
  G4MuonVDNuclearModel* mvdn = new G4MuonVDNuclearModel;


  G4ProcessManager* procMan = 0;

  // Gamma
  procMan = G4Gamma::Gamma()->GetProcessManager();
#if G4VERSION_NUMBER < 1100
  G4PhotoNuclearProcess* pnProc = new G4PhotoNuclearProcess;
#else
  auto* pnProc = new G4HadronInelasticProcess("PhotoNuclearProcess",
                             G4Gamma::Gamma() );
  pnProc->AddDataSet(new G4PhotoNuclearCrossSection);
#endif
  pnProc->RegisterMe(theGammaReaction);
  pnProc->RegisterMe(qgsp);
  procMan->AddDiscreteProcess(pnProc);

//#if G4VERSION_NUMBER >= 1100
//  auto* photonCapture = new G4HadronicProcess( "photonNuclear", fCapture );
//  auto* photonFission = new G4HadronicProcess( "photonFission", fFission );
//  procMan->AddDiscreteProcess(photonCapture);
//  procMan->AddDiscreteProcess(photonFission);
//#endif

  // Electron
  procMan = G4Electron::Electron()->GetProcessManager();
  G4ElectronNuclearProcess* emn = new G4ElectronNuclearProcess;
  emn->RegisterMe(evdn);
  procMan->AddDiscreteProcess(emn);

  // Positron
  procMan = G4Positron::Positron()->GetProcessManager();
  G4PositronNuclearProcess* epn = new G4PositronNuclearProcess;
  epn->RegisterMe(evdn);
  procMan->AddDiscreteProcess(epn);

  // Muon-
  procMan = G4MuonMinus::MuonMinus()->GetProcessManager();
  G4MuonNuclearProcess* mun = new G4MuonNuclearProcess;
  mun->RegisterMe(mvdn);
  procMan->AddDiscreteProcess(mun);

  // Muon+
  procMan = G4MuonPlus::MuonPlus()->GetProcessManager();
  procMan->AddDiscreteProcess(mun);
  
}


void eASTGammaLeptoNuclearPhysics::ConstructParticle()
{}

