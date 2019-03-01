// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        NeutronPhysics.cc                                            //
//  Description: Neutron hadronic physics constructor for EICPhysicsList      //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        3 July 2018                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "NeutronPhysics.hh"

#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4NeutronInelasticProcess.hh>
#include <Geant4/G4HadronElasticProcess.hh>
#include <Geant4/G4HadronCaptureProcess.hh>
#include <Geant4/G4NeutronKiller.hh>

#include <Geant4/G4CascadeInterface.hh>
#include <Geant4/G4TheoFSGenerator.hh>
#include <Geant4/G4FTFModel.hh>
#include <Geant4/G4ExcitedStringDecay.hh>
#include <Geant4/G4LundStringFragmentation.hh>
#include <Geant4/G4GeneratorPrecompoundInterface.hh>
#include <Geant4/G4ChipsElasticModel.hh>
#include <Geant4/G4NeutronRadCapture.hh>

#include <Geant4/G4BGGNucleonInelasticXS.hh>
#include <Geant4/G4NeutronElasticXS.hh>
#include <Geant4/G4NeutronCaptureXS.hh>

#include <Geant4/G4SystemOfUnits.hh>


NeutronPhysics::NeutronPhysics():
  ftfp(nullptr),
  stringModel(nullptr),
  stringDecay(nullptr),
  fragModel(nullptr),
  preCompoundModel(nullptr)
{}


NeutronPhysics::~NeutronPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}


void NeutronPhysics::ConstructParticle()
{}


void NeutronPhysics::ConstructProcess()
{
  // Low energy elastic model
  G4ChipsElasticModel* elMod = new G4ChipsElasticModel();

  // Use Bertini cascade for low energies
  G4CascadeInterface* loInelModel = new G4CascadeInterface;
  loInelModel->SetMinEnergy(0.0);
  loInelModel->SetMaxEnergy(12.0*GeV);

  // Capture model
  G4NeutronRadCapture* capModel = new G4NeutronRadCapture();

  // Use FTFP for high energies   ==>>   eventually replace this with new class FTFPInterface
  ftfp = new G4TheoFSGenerator("FTFP");
  stringModel = new G4FTFModel;
  stringDecay =
    new G4ExcitedStringDecay(fragModel = new G4LundStringFragmentation);
  stringModel->SetFragmentationModel(stringDecay);
  preCompoundModel = new G4GeneratorPrecompoundInterface();

  ftfp->SetHighEnergyGenerator(stringModel);
  ftfp->SetTransport(preCompoundModel); 
  ftfp->SetMinEnergy(5*GeV);
  ftfp->SetMaxEnergy(100*TeV);

  // Cross section sets
  G4BGGNucleonInelasticXS* inelCS = new G4BGGNucleonInelasticXS(G4Neutron::Neutron() );
  G4NeutronElasticXS* elCS = new G4NeutronElasticXS;
  G4NeutronCaptureXS* capCS = new G4NeutronCaptureXS;

  G4ProcessManager* procMan = G4Neutron::Neutron()->GetProcessManager();

  // Elastic process
  G4HadronElasticProcess* nProcEl = new G4HadronElasticProcess;
  nProcEl->RegisterMe(elMod);
  nProcEl->AddDataSet(elCS);
  procMan->AddDiscreteProcess(nProcEl);

  // Inelastic process
  G4NeutronInelasticProcess* nProcInel = new G4NeutronInelasticProcess;
  nProcInel->RegisterMe(loInelModel);
  nProcInel->RegisterMe(ftfp);
  nProcInel->AddDataSet(inelCS);
  procMan->AddDiscreteProcess(nProcInel);

  // Capture process
  G4HadronCaptureProcess* nProcCap = new G4HadronCaptureProcess("nCapture");
  nProcCap->RegisterMe(capModel);
  nProcCap->AddDataSet(capCS);
  procMan->AddDiscreteProcess(nProcCap);

  // Neutron cut (kill neutrons that live too long or have too little energy)
  G4NeutronKiller* nKiller = new G4NeutronKiller();
  nKiller->SetKinEnergyLimit(0.0*MeV);
  nKiller->SetTimeLimit(10.*microsecond);
  procMan->AddDiscreteProcess(nKiller);
  
}


