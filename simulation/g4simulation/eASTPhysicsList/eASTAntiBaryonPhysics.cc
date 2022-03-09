////////////////////////////////////////////////////////////////////////////////
//
//  eASTAntiBaryonPhysics.cc
//  Anti-baryon hadronic physics constructor for eASTPhysicsList
//
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.02.2021 : migration to Geant4 version 10.7 - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////


#include "eASTAntiBaryonPhysics.hh"

#include "G4ProcessManager.hh"

#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4AntiDeuteronInelasticProcess.hh"
#include "G4AntiTritonInelasticProcess.hh"
#include "G4AntiHe3InelasticProcess.hh"
#include "G4AntiAlphaInelasticProcess.hh"
#else
#include "G4HadronInelasticProcess.hh"
#endif

#include "G4HadronElasticProcess.hh"

#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4HadronElastic.hh"
#include "G4AntiNuclElastic.hh"
#include "G4HadronicAbsorptionFritiof.hh"

#include "G4ChipsAntiBaryonElasticXS.hh"
#include "G4ChipsHyperonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ChipsAntiBaryonElasticXS.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4CrossSectionElastic.hh"

#include "G4SystemOfUnits.hh"

#if G4VERSION_NUMBER < 1100
eASTAntiBaryonPhysics::eASTAntiBaryonPhysics()
{}

eASTAntiBaryonPhysics::~eASTAntiBaryonPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;

  delete theAntiNucleonXS;
}
#else
eASTAntiBaryonPhysics::eASTAntiBaryonPhysics()
: G4VPhysicsConstructor("eASTAntiBaryon")
{;}

eASTAntiBaryonPhysics::~eASTAntiBaryonPhysics()
{;}
#endif

void eASTAntiBaryonPhysics::ConstructParticle()
{}


void eASTAntiBaryonPhysics::ConstructProcess()
{
  G4ProcessManager* procMan = 0;

  // One elastic model for all anti-hyperon and anti-neutron energies
  G4HadronElastic* elModel = new G4HadronElastic();

  // Elastic models for anti-(p, d, t, He3, alpha)  
  G4HadronElastic* loelModel = new G4HadronElastic();
  loelModel->SetMaxEnergy(100.1*MeV);

  G4AntiNuclElastic* anucEl = new G4AntiNuclElastic();
  anucEl->SetMinEnergy(100.0*MeV);

  // Use FTFP for all energies   ==>>   eventually replace this with new class FTFPInterface
  ftfp = new G4TheoFSGenerator("FTFP");
  stringModel = new G4FTFModel;
  stringDecay =
    new G4ExcitedStringDecay(fragModel = new G4LundStringFragmentation);
  stringModel->SetFragmentationModel(stringDecay);
  preCompoundModel = new G4GeneratorPrecompoundInterface();

  ftfp->SetHighEnergyGenerator(stringModel);
  ftfp->SetTransport(preCompoundModel); 
  ftfp->SetMinEnergy(0.0);
  ftfp->SetMaxEnergy(100*TeV);

  // Elastic data sets
  G4CrossSectionElastic* anucElxs =
    new G4CrossSectionElastic(anucEl->GetComponentCrossSection() );
  G4VCrossSectionDataSet* abaryElXs = new G4ChipsAntiBaryonElasticXS;

  G4VCrossSectionDataSet* anucnucElxs =
    new G4CrossSectionElastic(new G4ComponentAntiNuclNuclearXS);

  // Inelastic cross section sets
  theAntiNucleonXS = new G4ComponentAntiNuclNuclearXS;
  G4VCrossSectionDataSet* antiNucleonData =
    new G4CrossSectionInelastic(theAntiNucleonXS);

  G4ChipsHyperonInelasticXS* hchipsInelastic = new G4ChipsHyperonInelasticXS;

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-proton                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiProton::AntiProton()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* apProcEl = new G4HadronElasticProcess;
  apProcEl->RegisterMe(loelModel);
  apProcEl->RegisterMe(anucEl);
  apProcEl->AddDataSet(anucElxs);
  procMan->AddDiscreteProcess(apProcEl);

  // inelastic 
#if G4VERSION_NUMBER < 1100
  G4AntiProtonInelasticProcess* apProcInel = new G4AntiProtonInelasticProcess;
#else
  auto* apProcInel = new G4HadronInelasticProcess("AntiProtonInelasticProcess",
                                 G4AntiProton::AntiProton() );
#endif
  apProcInel->RegisterMe(ftfp);
  apProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(apProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* apAbsorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(apAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-neutron                                                           // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* anProcEl = new G4HadronElasticProcess;
  anProcEl->RegisterMe(elModel);
  anProcEl->AddDataSet(anucnucElxs);
  procMan->AddDiscreteProcess(anProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiNeutronInelasticProcess* anProcInel = new G4AntiNeutronInelasticProcess;
#else
  auto* anProcInel = new G4HadronInelasticProcess("AntiNeutronInelasticProcess",
                                 G4AntiNeutron::AntiNeutron() );
#endif
  anProcInel->RegisterMe(ftfp);
  anProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(anProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-deuteron                                                          // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiDeuteron::AntiDeuteron()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* adProcEl = new G4HadronElasticProcess;
  adProcEl->RegisterMe(loelModel);
  adProcEl->RegisterMe(anucEl);
  adProcEl->AddDataSet(anucElxs);
  procMan->AddDiscreteProcess(adProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiDeuteronInelasticProcess* adProcInel = new G4AntiDeuteronInelasticProcess;
#else
  auto* adProcInel = new G4HadronInelasticProcess("AntiDeuteronInelasticProcess",
                                 G4AntiDeuteron::AntiDeuteron() );
#endif
  adProcInel->RegisterMe(ftfp);
  adProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(adProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* adAbsorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(adAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-triton                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiTriton::AntiTriton()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* atProcEl = new G4HadronElasticProcess;
  atProcEl->RegisterMe(loelModel);
  atProcEl->RegisterMe(anucEl);
  atProcEl->AddDataSet(anucElxs);
  procMan->AddDiscreteProcess(atProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiTritonInelasticProcess* atProcInel = new G4AntiTritonInelasticProcess;
#else
  auto* atProcInel = new G4HadronInelasticProcess("AntiTritonInelasticProcess",
                                 G4AntiTriton::AntiTriton() );
#endif
  atProcInel->RegisterMe(ftfp);
  atProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(atProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* atAbsorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(atAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-He3                                                               // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiHe3::AntiHe3()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* ahe3ProcEl = new G4HadronElasticProcess;
  ahe3ProcEl->RegisterMe(loelModel);
  ahe3ProcEl->RegisterMe(anucEl);
  ahe3ProcEl->AddDataSet(anucElxs);
  procMan->AddDiscreteProcess(ahe3ProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiHe3InelasticProcess* ahe3ProcInel = new G4AntiHe3InelasticProcess;
#else
  auto* ahe3ProcInel = new G4HadronInelasticProcess("Anti3HeInelasticProcess",
                                   G4AntiHe3::AntiHe3() );
#endif
  ahe3ProcInel->RegisterMe(ftfp);
  ahe3ProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(ahe3ProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* ahe3Absorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(ahe3Absorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-alpha                                                             // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiAlpha::AntiAlpha()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* aaProcEl = new G4HadronElasticProcess;
  aaProcEl->RegisterMe(loelModel);
  aaProcEl->RegisterMe(anucEl);
  aaProcEl->AddDataSet(anucElxs);
  procMan->AddDiscreteProcess(aaProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiAlphaInelasticProcess* aaProcInel = new G4AntiAlphaInelasticProcess;
#else
  auto* aaProcInel = new G4HadronInelasticProcess("AntiAlphaInelasticProcess",
                                 G4AntiAlpha::AntiAlpha() );
#endif
  aaProcInel->RegisterMe(ftfp);
  aaProcInel->AddDataSet(antiNucleonData);
  procMan->AddDiscreteProcess(aaProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* aaAbsorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(aaAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-lambda                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiLambda::AntiLambda()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* alamProcEl = new G4HadronElasticProcess;
  alamProcEl->RegisterMe(elModel);
  alamProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(alamProcEl);

  // inelastic 
#if G4VERSION_NUMBER < 1100
  G4AntiLambdaInelasticProcess* alamProcInel = new G4AntiLambdaInelasticProcess;
#else
  auto* alamProcInel = new G4HadronInelasticProcess("AntiLambdaInelasticProcess",
                                   G4AntiLambda::AntiLambda() );
#endif
  alamProcInel->RegisterMe(ftfp);
  alamProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(alamProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-sigma+                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* aspProcEl = new G4HadronElasticProcess;
  aspProcEl->RegisterMe(elModel);
  aspProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(aspProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiSigmaPlusInelasticProcess* aspProcInel = new G4AntiSigmaPlusInelasticProcess;
#else
  auto* aspProcInel = new G4HadronInelasticProcess("AntiSigmaPInelasticProcess",
                                  G4AntiSigmaPlus::AntiSigmaPlus() );
#endif
  aspProcInel->RegisterMe(ftfp);
  aspProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(aspProcInel);

  // stopping
  G4HadronicAbsorptionFritiof* aspAbsorb = new G4HadronicAbsorptionFritiof();
  procMan->AddRestProcess(aspAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-sigma-                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* asmProcEl = new G4HadronElasticProcess;
  asmProcEl->RegisterMe(elModel);
  asmProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(asmProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiSigmaMinusInelasticProcess* asmProcInel = new G4AntiSigmaMinusInelasticProcess;
#else
  auto* asmProcInel = new G4HadronInelasticProcess("AntiSigmaMInelasticProcess",
                                  G4AntiSigmaMinus::AntiSigmaMinus() );
#endif
  asmProcInel->RegisterMe(ftfp);
  asmProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(asmProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-xi0                                                               // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* axzProcEl = new G4HadronElasticProcess;
  axzProcEl->RegisterMe(elModel);
  axzProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(axzProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiXiZeroInelasticProcess* axzProcInel = new G4AntiXiZeroInelasticProcess;
#else
  auto* axzProcInel = new G4HadronInelasticProcess("AntiXi0InelasticProcess",
                                  G4AntiXiZero::AntiXiZero() );
#endif
  axzProcInel->RegisterMe(ftfp);
  axzProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(axzProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-xi-                                                               // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* axmProcEl = new G4HadronElasticProcess;
  axmProcEl->RegisterMe(elModel);
  axmProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(axmProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiXiMinusInelasticProcess* axmProcInel = new G4AntiXiMinusInelasticProcess;
#else
  auto* axmProcInel = new G4HadronInelasticProcess("AntiXiMInelasticProcess",
                                  G4AntiXiMinus::AntiXiMinus() );
#endif
  axmProcInel->RegisterMe(ftfp);
  axmProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(axmProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   Anti-omega-                                                            // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* aomProcEl = new G4HadronElasticProcess;
  aomProcEl->RegisterMe(elModel);
  aomProcEl->AddDataSet(abaryElXs);
  procMan->AddDiscreteProcess(aomProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4AntiOmegaMinusInelasticProcess* aomProcInel = new G4AntiOmegaMinusInelasticProcess;
#else
  auto* aomProcInel = new G4HadronInelasticProcess("AntiOmegaMInelasticProcess",
                                  G4AntiOmegaMinus::AntiOmegaMinus() );
#endif
  aomProcInel->RegisterMe(ftfp);
  aomProcInel->AddDataSet(hchipsInelastic);
  procMan->AddDiscreteProcess(aomProcInel);

}

void eASTAntiBaryonPhysics::TerminateWorker()
{}

