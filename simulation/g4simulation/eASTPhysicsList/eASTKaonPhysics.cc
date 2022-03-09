////////////////////////////////////////////////////////////////////////////////
//
//  eASTKaonPhysics.hh
//  Kaon hadronic physics constructor for eASTPhysicsList
//
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC)
//    May.02.2021 : migration to Genat4 version 10.7 - Dennis H. Wright (SLAC)
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//    Dec.22.2021 : migration to Geant4 version 11.0 - Makoto Asai (JLab)
//
////////////////////////////////////////////////////////////////////////////////


#include "eASTKaonPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1100
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
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

#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionElastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4SystemOfUnits.hh"

#if G4VERSION_NUMBER < 1100
eASTKaonPhysics::eASTKaonPhysics()
{}

eASTKaonPhysics::~eASTKaonPhysics()
{
  delete stringDecay;
  delete stringModel;
  delete fragModel;
  delete preCompoundModel;
}
#else
eASTKaonPhysics::eASTKaonPhysics()
: G4VPhysicsConstructor("eASTKaon")
{;}

eASTKaonPhysics::~eASTKaonPhysics()
{;}
#endif

void eASTKaonPhysics::ConstructParticle()
{}


void eASTKaonPhysics::ConstructProcess()
{
  G4ProcessManager* procMan;

  // One elastic model for all kaon energies
  G4HadronElastic* elModel = new G4HadronElastic();

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

  // Inelastic cross section sets
  G4VCrossSectionDataSet* kpCS = new G4ChipsKaonPlusInelasticXS;
  G4VCrossSectionDataSet* kmCS = new G4ChipsKaonMinusInelasticXS;
  G4VCrossSectionDataSet* kzCS = new G4ChipsKaonZeroInelasticXS;

  // Elastic cross section
  G4VCrossSectionDataSet* kelCS =
    new G4CrossSectionElastic(new G4ComponentGGHadronNucleusXsc);

  //////////////////////////////////////////////////////////////////////////////
  //   K+                                                                     // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4KaonPlus::KaonPlus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* kpProcEl = new G4HadronElasticProcess;
  kpProcEl->RegisterMe(elModel);
  kpProcEl->AddDataSet(kelCS);
  procMan->AddDiscreteProcess(kpProcEl);

  // inelastic 
#if G4VERSION_NUMBER < 1100
  G4KaonPlusInelasticProcess* kpProcInel = new G4KaonPlusInelasticProcess;
#else 
  auto* kpProcInel = new G4HadronInelasticProcess("KaonPlusInelasticProcess",
                                 G4KaonPlus::KaonPlus() );
#endif
  kpProcInel->RegisterMe(loInelModel);
  kpProcInel->RegisterMe(ftfp);
  kpProcInel->AddDataSet(kpCS);
  procMan->AddDiscreteProcess(kpProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   K-                                                                     // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4KaonMinus::KaonMinus()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* kmProcEl = new G4HadronElasticProcess;
  kmProcEl->RegisterMe(elModel);
  kmProcEl->AddDataSet(kelCS);
  procMan->AddDiscreteProcess(kmProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4KaonMinusInelasticProcess* kmProcInel = new G4KaonMinusInelasticProcess;
#else
  auto* kmProcInel = new G4HadronInelasticProcess("KaonMinusInelasticProcess",
                                 G4KaonMinus::KaonMinus() );
#endif
  kmProcInel->RegisterMe(loInelModel);
  kmProcInel->RegisterMe(ftfp);
  kmProcInel->AddDataSet(kmCS);
  procMan->AddDiscreteProcess(kmProcInel);

  // stopping
  G4HadronicAbsorptionBertini* bertAbsorb = new G4HadronicAbsorptionBertini;
  procMan->AddRestProcess(bertAbsorb);

  //////////////////////////////////////////////////////////////////////////////
  //   K0L                                                                    // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* k0LProcEl = new G4HadronElasticProcess;
  k0LProcEl->RegisterMe(elModel);
  k0LProcEl->AddDataSet(kelCS);
  procMan->AddDiscreteProcess(k0LProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4KaonZeroLInelasticProcess* k0LProcInel = new G4KaonZeroLInelasticProcess;
#else
  auto* k0LProcInel = new G4HadronInelasticProcess("Kaon0LongInelasticProcess",
                                  G4KaonZeroLong::KaonZeroLong() );
#endif
  k0LProcInel->RegisterMe(loInelModel);
  k0LProcInel->RegisterMe(ftfp);
  k0LProcInel->AddDataSet(kzCS);
  procMan->AddDiscreteProcess(k0LProcInel);

  //////////////////////////////////////////////////////////////////////////////
  //   K0S                                                                    // 
  //////////////////////////////////////////////////////////////////////////////

  procMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();

  // elastic
  G4HadronElasticProcess* k0SProcEl = new G4HadronElasticProcess;
  k0SProcEl->RegisterMe(elModel);
  k0SProcEl->AddDataSet(kelCS);
  procMan->AddDiscreteProcess(k0SProcEl);

  // inelastic
#if G4VERSION_NUMBER < 1100
  G4KaonZeroSInelasticProcess* k0SProcInel = new G4KaonZeroSInelasticProcess;
#else
  auto* k0SProcInel = new G4HadronInelasticProcess("Kaon0ShortInelasticProcess",
                                  G4KaonZeroShort::KaonZeroShort() );
#endif
  k0SProcInel->RegisterMe(loInelModel);
  k0SProcInel->RegisterMe(ftfp);
  k0SProcInel->AddDataSet(kzCS);
  procMan->AddDiscreteProcess(k0SProcInel);
}

void eASTKaonPhysics::TerminateWorker()
{}

