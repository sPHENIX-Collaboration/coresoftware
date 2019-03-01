// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        GammaLeptoNuclearPhysics.cc                                  //
//  Description: Gamma-nuclear, electro-nuclear and muon-nuclear physics      // 
//                 constructor for EICPhysicsList                             //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        20 July 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "GammaLeptoNuclearPhysics.hh"

#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4PhotoNuclearProcess.hh>
#include <Geant4/G4ElectronNuclearProcess.hh>
#include <Geant4/G4PositronNuclearProcess.hh>
#include <Geant4/G4MuonNuclearProcess.hh>

#include <Geant4/G4CascadeInterface.hh>
#include <Geant4/G4ElectroVDNuclearModel.hh>
#include <Geant4/G4MuonVDNuclearModel.hh>

#include <Geant4/G4TheoFSGenerator.hh>
#include <Geant4/G4ExcitedStringDecay.hh>
#include <Geant4/G4QGSMFragmentation.hh>
#include <Geant4/G4GeneratorPrecompoundInterface.hh>

#include <Geant4/G4SystemOfUnits.hh>


GammaLeptoNuclearPhysics::GammaLeptoNuclearPhysics():
  qgsp(nullptr),
  stringModel(nullptr),
  stringDecay(nullptr),
  fragModel(nullptr),
  preCompoundModel(nullptr)
{}


GammaLeptoNuclearPhysics::~GammaLeptoNuclearPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}


void GammaLeptoNuclearPhysics::ConstructProcess()
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
  G4PhotoNuclearProcess* pnProc = new G4PhotoNuclearProcess;
  pnProc->RegisterMe(theGammaReaction);
  pnProc->RegisterMe(qgsp);
  procMan->AddDiscreteProcess(pnProc);

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


void GammaLeptoNuclearPhysics::ConstructParticle()
{}

