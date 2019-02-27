// $Id: $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        ProtonPhysics.cc                                             //
//  Description: Proton hadronic physics constructor for EICPhysicsList       //
//                                                                            //
//  Author:      Dennis H. Wright (SLAC)                                      //  
//  Date:        22 June 2018                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "ProtonPhysics.hh"

#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4ProtonInelasticProcess.hh>
#include <Geant4/G4HadronElasticProcess.hh>

#include <Geant4/G4CascadeInterface.hh>
#include <Geant4/G4TheoFSGenerator.hh>
#include <Geant4/G4FTFModel.hh>
#include <Geant4/G4ExcitedStringDecay.hh>
#include <Geant4/G4LundStringFragmentation.hh>
#include <Geant4/G4GeneratorPrecompoundInterface.hh>
#include <Geant4/G4ChipsElasticModel.hh>

#include <Geant4/G4BGGNucleonInelasticXS.hh>
#include <Geant4/G4ChipsProtonElasticXS.hh>

#include <Geant4/G4SystemOfUnits.hh>


ProtonPhysics::ProtonPhysics()
{}


ProtonPhysics::~ProtonPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}


void ProtonPhysics::ConstructProcess()
{
  // Low energy elastic model
  G4ChipsElasticModel* elMod = new G4ChipsElasticModel();

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
  ftfp->SetMinEnergy(5*GeV);
  ftfp->SetMaxEnergy(100*TeV);

  // Inelastic cross section
  G4BGGNucleonInelasticXS* inelCS = new G4BGGNucleonInelasticXS(G4Proton::Proton() );
  G4ChipsProtonElasticXS* elCS = new G4ChipsProtonElasticXS;

  G4ProcessManager* procMan = G4Proton::Proton()->GetProcessManager();

  // Elastic process
  G4HadronElasticProcess* pProcEl = new G4HadronElasticProcess;
  pProcEl->RegisterMe(elMod);
  pProcEl->AddDataSet(elCS);
  procMan->AddDiscreteProcess(pProcEl);

  // Inelastic process
  G4ProtonInelasticProcess* pProcInel = new G4ProtonInelasticProcess;
  pProcInel->RegisterMe(loInelModel);
  pProcInel->RegisterMe(ftfp);
  pProcInel->AddDataSet(inelCS);
  procMan->AddDiscreteProcess(pProcInel);

}


void ProtonPhysics::ConstructParticle()
{}

