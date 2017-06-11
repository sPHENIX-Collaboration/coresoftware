#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

#include <TPCbase/TPCHit.h>
#include <TPCbase/TPCHitsContainer.h>
#include <TPCbase/TPCDigitsContainer.h>
#include <TPCbase/TPCPadMap.h>
#include <TPCbase/TPCDataTypes.h>
#include <TPCbase/TPCConstants.h>
#include "TPCSimulation.h"
#include "TPCCloud.h"

using namespace TPCDataTypes;
using namespace TPCConstants;

//=====
TPCSimulation::TPCSimulation() :
  fHits(NULL),
  fDigits(NULL),
  fPadMap( new TPCPadMap ),
  fCloud( new TPCCloud ),
  fVerbosity(0),
  fTreeFileName(""),
  fTree(NULL),
  fFile(NULL)
{
}
//=====
TPCSimulation::~TPCSimulation() {
  delete fPadMap;
  delete fCloud;
  if(fFile) delete fFile;
}
//====
void TPCSimulation::Hits2Digits() {
  if(!fHits) return;
  if(!fDigits) return;
  Hits2Digits(0);
  Hits2Digits(1);
  return;
}
//=====
void TPCSimulation::Edep2Ele()
{
  // converts ionization energy
  // into total electrons
  if(fVerbosity>2) std::cout << "Edep2Ele ... begin." << std::endl;
  Float_t ene = fCloud->GetDEnergy();
  ene = gRandom->PoissonD( ene*kElectronsPerKeV );
  fCloud->SetDEnergy(ene);
  if(fVerbosity>2) std::cout << "Edep2Ele ... end." << std::endl;
}
//=====
void TPCSimulation::Transport()
{
  // compute path distortion due
  // to fields and ion backflow
  if(fVerbosity>2) std::cout << "Transport ... begin." << std::endl;
  Float_t dr = 0.0;
  Float_t dz = 0.0;
  Float_t dp = 0.0;
  fCloud->AddR( dr );    // radial distortion
  fCloud->AddPhi( dp );  // phi distortion
  fCloud->AddZ( dz );    // z distortion
  fCloud->AddMST( 0.0 ); // add traverse variance
  fCloud->AddMSL0( 0.0 ); // add longitudinal variance
  // compute diffusion spread
  Float_t dlength = kGasHalfOfLength - TMath::Abs(fCloud->GetZ()+dz); // TODO: add path distortion in r and phi
  Float_t tt = kGasDiffusionTransverseSquare*dlength;
  //Float_t ll = kGasDiffusionLongitudinalSquare*dlength;
  fCloud->AddMST( tt );
  //fCloud->AddMSL0( ll );
  //fCloud->AddMSL1( ll );
  if(fVerbosity>2) std::cout << "Transport ... end." << std::endl;
}
//====
void TPCSimulation::Amplify() {
  fCloud->SetElectrons(2000*fCloud->GetDEnergy());
  fCloud->AddMST( 0.01 ); // GEM smeering
}
//====
void TPCSimulation::Digitize() {
  Float_t dlength = kGasHalfOfLength-TMath::Abs(fCloud->GetZ());
  fCloud->SetZ( dlength*kInverseOfDriftVelocity ); // moves z(cm) to time
  fCloud->SetMSL0( kTimeShapeRiseSquare );
  fCloud->SetMSL1( kTimeShapeTailSquare );
}
//====
void TPCSimulation::Hits2Digits(int blk) {
  UInt_t nhits = fHits->GetNHits(blk);
  if(fVerbosity>1) std::cout << "TPCSimulation::Hits2Digits(" << blk << ") => hits: " << nhits << std::endl;
  TPCHit *hit;
  for(UInt_t i=0; i!=nhits; ++i) {
    if(fVerbosity>1)
      if(int(i*100.0/nhits)%10==0) std::cout << i << "/" << nhits << std::endl;
    hit = fHits->GetHit(blk,i);
    if(!hit) continue;
    if(fTree) {
      fHitRegister.track = hit->GetTrack();
      fHitRegister.r = hit->GetR();
      fHitRegister.phi = hit->GetPhi();
      fHitRegister.z = hit->GetZ();
      fHitRegister.length = hit->GetL();
      fHitRegister.edep = hit->GetDEnergy();
      fTree->Fill();
    }
    fCloud->CopyFrom( hit );
    Edep2Ele();
    Transport();
    Amplify();
    Digitize();
    if(fVerbosity>3) std::cout << "rmst " << fCloud->GetRMST() << " rmsl0 " << fCloud->GetRMSL0()
			       << " rmsl1 " << fCloud->GetRMSL1() << std::endl;
    ModuleRange_t mods = fPadMap->FindModuleRP( fCloud->GetR(), fCloud->GetPhi(),
						fCloud->GetRMST(), fCloud->GetElectrons() );
    if(fVerbosity>2) std::cout << "range of modules " << mods.size() << std::endl;
    for(unsigned int j=0; j!=mods.size(); ++j)
      PushCloud2Module(mods[j]+blk*kNModulesPerPlate);
  }
  if(fVerbosity>1) std::cout << "TPCSimulation::Hits2Digits(" << blk << ") => pushed " << std::endl;
}
//=====
void TPCSimulation::PushCloud2Module(Module_t mod)
{
  if(fVerbosity>1) std::cout << "PushCloud2Module " << int(mod) << std::endl;
  Float_t rmst = fCloud->GetRMST();
  PadQuotaRange_t qRange = fPadMap->GetPadQuotasRP( fCloud->GetElectrons(), mod, fCloud->GetR(), fCloud->GetPhi(), rmst );
  if(qRange.size()<1) return; // no coincidence found
  Float_t rmsl0 = fCloud->GetRMSL0();
  Float_t rmsl1 = fCloud->GetRMSL1();
  Pad_t prc;
  Time_t tim;
  Float_t iww;
  Adc_t jww;
  if(fVerbosity>3) {
  std::cout << "CLOUD R " << fCloud->GetR() << " PHI " << fCloud->GetPhi() << " T " << fCloud->GetZ() << " || "
	    << " RMST " << rmst << " RMSL " << rmsl0 << " " << rmsl1 << std::endl;
  }
  if(fVerbosity>1) std::cout << "CLOUD WEIGHT " << fCloud->GetElectrons() <<
		     " to be stored in " << qRange.size() <<  " pads." << std::endl;
  Float_t sumiww = 0;
  for(unsigned int i=0; i!=qRange.size(); ++i) { // loop over traverse range
    prc = qRange[i].first;
    iww = qRange[i].second;
    sumiww += iww;
    TimeQuotaRange_t tRange = fPadMap->GetTimeQuotas( iww, fCloud->GetZ(), rmsl0, rmsl1 );
    if(fVerbosity>1) std::cout << "  PAD " << int(prc) << "(" << fPadMap->GetRow(mod,prc) << ")" << " weight "
			       << iww << " to be stored in " << tRange.size() << " time bins." << std::endl;
    if(fVerbosity>3) {
      PairOfFloats_t rp = fPadMap->GetRP(mod,prc);
      PairOfFloats_t dxy = fPadMap->GetDXDY(mod,prc);
      std::cout << "  R " << rp.first << " PHI " << rp.second
		<< " DX "  << dxy.first << " DY " << dxy.second << std::endl;
    }
    Adc_t sumjww = 0;
    for(unsigned int j=0; j!=tRange.size(); ++j) { // loop over time axis
      tim = tRange[j].first;
      jww = tRange[j].second;
      sumjww += jww;
      if(fVerbosity>2) std::cout << "    TIME " << int(tim) << " weight " << jww << std::endl;
      fDigits->Add( mod, prc, tim, jww );
    }
    if(fVerbosity>2) std::cout << "    ==> SUM TIME WIGHTS " << int(sumjww) << std::endl;
  }
  if(fVerbosity>1) std::cout << "  ==> SUM PAD WIGHTS " << int(sumiww) << std::endl;
}
//=====
void TPCSimulation::PrepareTree(TString name)
{
  fFile = new TFile(name.Data(),"RECREATE");
  fTree = new TTree("TPCSimData","TPCSimData");
  fTree->SetAutoSave(10000);
  fTree->Branch("hits",&fHitRegister,"track/I:r/F:phi:z:length:edep");
}
//=====
void TPCSimulation::WriteFile()
{
  std::cout << "TPCSimulation::Writing into file : " << fFile->GetNewUrl() << std::endl;
  fFile->cd();
  fFile->Write();
  fFile->Close();
}
