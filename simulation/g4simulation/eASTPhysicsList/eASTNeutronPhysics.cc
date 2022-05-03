////////////////////////////////////////////////////////////////////////////////
//
//  eASTNeutronPhysics.cc
//  Neutron hadronic physics constructor for eASTPhysicsList
//
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////


#include "eASTNeutronPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#else
#include "G4HadronInelasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#endif
#include "G4HadronElasticProcess.hh"
#include "G4NeutronKiller.hh"

#include "G4CascadeInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ChipsElasticModel.hh"
#include "G4NeutronRadCapture.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4NeutronElasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4SystemOfUnits.hh"

#if G4VERSION_NUMBER < 1100
eASTNeutronPhysics::eASTNeutronPhysics()
{
}

eASTNeutronPhysics::~eASTNeutronPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}
#else
eASTNeutronPhysics::eASTNeutronPhysics()
: G4VPhysicsConstructor("eASTNeutron")
{;}

eASTNeutronPhysics::~eASTNeutronPhysics()
{;}
#endif

void eASTNeutronPhysics::ConstructParticle()
{}


void eASTNeutronPhysics::ConstructProcess()
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
#if G4VERSION_NUMBER < 1100
  G4NeutronInelasticProcess* nProcInel = new G4NeutronInelasticProcess;
#else
  auto* nProcInel = new G4HadronInelasticProcess("NeutronInelasticProcess",
                                G4Neutron::Neutron() );
#endif
  nProcInel->RegisterMe(loInelModel);
  nProcInel->RegisterMe(ftfp);
  nProcInel->AddDataSet(inelCS);
  procMan->AddDiscreteProcess(nProcInel);

  // Capture process
#if G4VERSION_NUMBER < 1100
  G4HadronCaptureProcess* nProcCap = new G4HadronCaptureProcess("nCapture");
#else
  auto* nProcCap = new G4NeutronCaptureProcess("nCapture");
#endif
  nProcCap->RegisterMe(capModel);
  nProcCap->AddDataSet(capCS);
  procMan->AddDiscreteProcess(nProcCap);

  // Neutron cut (kill neutrons that live too long or have too little energy)
  G4NeutronKiller* nKiller = new G4NeutronKiller();
  nKiller->SetKinEnergyLimit(0.0*MeV);
  nKiller->SetTimeLimit(10.*microsecond);
  procMan->AddDiscreteProcess(nKiller);
  
}


