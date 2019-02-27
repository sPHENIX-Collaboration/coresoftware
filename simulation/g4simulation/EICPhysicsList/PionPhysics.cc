// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        PionPhysics.cc                                               //
//  Description: Pion hadronic physics constructor for EICPhysicsList         //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        21 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "PionPhysics.hh"

#include <Geant4/G4MesonConstructor.hh>

#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4PionPlusInelasticProcess.hh>
#include <Geant4/G4PionMinusInelasticProcess.hh>
#include <Geant4/G4HadronElasticProcess.hh>
#include <Geant4/G4HadronicAbsorptionBertini.hh>

#include <Geant4/G4CascadeInterface.hh>
#include <Geant4/G4TheoFSGenerator.hh>
#include <Geant4/G4FTFModel.hh>
#include <Geant4/G4ExcitedStringDecay.hh>
#include <Geant4/G4LundStringFragmentation.hh>
#include <Geant4/G4GeneratorPrecompoundInterface.hh>
#include <Geant4/G4HadronElastic.hh>
#include <Geant4/G4ElasticHadrNucleusHE.hh>

#include <Geant4/G4PiNuclearCrossSection.hh>
#include <Geant4/G4CrossSectionPairGG.hh>
#include <Geant4/G4BGGPionElasticXS.hh>

#include <Geant4/G4SystemOfUnits.hh>


PionPhysics::PionPhysics()
{}


PionPhysics::~PionPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}


void PionPhysics::ConstructParticle()
{}


void PionPhysics::ConstructProcess()
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

  // Inelastic cross section
  G4VCrossSectionDataSet* piCS = new G4CrossSectionPairGG(new G4PiNuclearCrossSection, 91*GeV);

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
  G4PionPlusInelasticProcess* pipProcInel = new G4PionPlusInelasticProcess;
  pipProcInel->RegisterMe(loInelModel);
  pipProcInel->RegisterMe(ftfp);
  pipProcInel->AddDataSet(piCS);
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
  G4PionMinusInelasticProcess* pimProcInel = new G4PionMinusInelasticProcess;
  pimProcInel->RegisterMe(loInelModel);
  pimProcInel->RegisterMe(ftfp);
  pimProcInel->AddDataSet(piCS);
  procMan->AddDiscreteProcess(pimProcInel);

  // stopping
  G4HadronicAbsorptionBertini* bertAbsorb = new G4HadronicAbsorptionBertini;
  procMan->AddRestProcess(bertAbsorb);

}

