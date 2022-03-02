////////////////////////////////////////////////////////////////////////////////
//
//  eASTPionPhysics.hh
//  Pion hadronic physics constructor for eASTPhysicsList
//
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.02.2021 : migration to Geant4 version 10.7 - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////


#include "eASTPionPhysics.hh"
#include "G4MesonConstructor.hh"

#include "G4ProcessManager.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#else
#include "G4HadronInelasticProcess.hh"
#endif
#include "G4HadronElasticProcess.hh"
#include "G4HadronicAbsorptionBertini.hh"

#include "G4CascadeInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4HadronElastic.hh"
#include "G4ElasticHadrNucleusHE.hh"

#include "G4BGGPionElasticXS.hh"
#include "G4BGGPionInelasticXS.hh"

#include "G4SystemOfUnits.hh"

#if G4VERSION_NUMBER < 1100
eASTPionPhysics::eASTPionPhysics()
{}

eASTPionPhysics::~eASTPionPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}
#else
eASTPionPhysics::eASTPionPhysics()
: G4VPhysicsConstructor("eASTPion")
{;}

eASTPionPhysics::~eASTPionPhysics()
{;}
#endif

void eASTPionPhysics::ConstructParticle()
{}


void eASTPionPhysics::ConstructProcess()
{
  G4ProcessManager* procMan;

  // Low energy elastic model
  G4HadronElastic* loElModel = new G4HadronElastic();
  loElModel->SetMaxEnergy(1.0001*GeV);

  // High energy elastic model
  G4ElasticHadrNucleusHE* hiElModel = new G4ElasticHadrNucleusHE();
  hiElModel->SetMinEnergy(1.0*GeV);

  // Use Bertini cascade for low energies
  G4CascadeInterface* loInelModel = new G4CascadeInterface;
  loInelModel->SetMinEnergy(0.0);
  loInelModel->SetMaxEnergy(12.0*GeV);

  // Use FTFP for high energies   ==>>   eventually replace this with new class FTFPInterface
  ftfp = new G4TheoFSGenerator("FTFP");
  stringModel = new G4FTFModel;
  stringDecay =
    new G4ExcitedStringDecay(fragModel = new G4LundStringFragmentation);
  stringModel->SetFragmentationModel(stringDecay);
  preCompoundModel = new G4GeneratorPrecompoundInterface();

  ftfp->SetHighEnergyGenerator(stringModel);
  ftfp->SetTransport(preCompoundModel); 
  ftfp->SetMinEnergy(10*GeV);
  ftfp->SetMaxEnergy(100*TeV);

  //////////////////////////////////////////////////////////////////////////////
  //   pi+                                                                    // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4PionPlus::PionPlus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* pipProcEl = new G4HadronElasticProcess;
  pipProcEl->RegisterMe(loElModel);
  pipProcEl->RegisterMe(hiElModel);
  pipProcEl->AddDataSet(new G4BGGPionElasticXS(G4PionPlus::PionPlus() ) );
  procMan->AddDiscreteProcess(pipProcEl);

  // inelastic 
#if G4VERSION_NUMBER < 1100
  G4PionPlusInelasticProcess* pipProcInel = new G4PionPlusInelasticProcess;
#else
  auto* pipProcInel = new G4HadronInelasticProcess("PionPlusInelasticProcess",
                                  G4PionPlus::PionPlus() );
#endif
  pipProcInel->RegisterMe(loInelModel);
  pipProcInel->RegisterMe(ftfp);
  pipProcInel->AddDataSet(new G4BGGPionInelasticXS(G4PionPlus::PionPlus() ) );
  procMan->AddDiscreteProcess(pipProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   pi-                                                                    // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4PionMinus::PionMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* pimProcEl = new G4HadronElasticProcess;
  pimProcEl->RegisterMe(loElModel);
  pimProcEl->RegisterMe(hiElModel);
  pimProcEl->AddDataSet(new G4BGGPionElasticXS(G4PionMinus::PionMinus() ) );
  procMan->AddDiscreteProcess(pimProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4PionMinusInelasticProcess* pimProcInel = new G4PionMinusInelasticProcess;
#else
  auto* pimProcInel = new G4HadronInelasticProcess("PionMinusInelasticProcess",
                                  G4PionMinus::PionMinus() );
#endif
  pimProcInel->RegisterMe(loInelModel);
  pimProcInel->RegisterMe(ftfp);
  pimProcInel->AddDataSet(new G4BGGPionInelasticXS(G4PionMinus::PionMinus() ) );
  procMan->AddDiscreteProcess(pimProcInel);

  // stopping
  G4HadronicAbsorptionBertini* bertAbsorb = new G4HadronicAbsorptionBertini;
  procMan->AddRestProcess(bertAbsorb);

}

