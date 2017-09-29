#include <TFile.h>
#include <TTree.h>

#ifndef __CINT__
#include <gsl/gsl_randist.h>
#include <phool/PHRandomSeed.h>
#else
#include <TRandom3.h>
#endif

#include <TPCHit.h>
#include <TPCHitsContainer.h>
#include <TPCDigitsContainer.h>
#include <TPCPadMap.h>
#include <TPCDataTypes.h>
#include <TPCConstants.h>
#include "TPCSimulation.h"
#include "TPCElectronPDF.h"

using namespace TPCDataTypes;
using namespace TPCConstants;

//=====
TPCSimulation::TPCSimulation() :
  fHits(NULL),
  fDigits(NULL),
  fPadMap( new TPCPadMap ),
  fEPDF( new TPCElectronPDF ),
  fVerbosity(0),
  fTreeFileName(""),
  fTree(NULL),
  fFile(NULL)
{
#ifndef __CINT__
  fRandom = gsl_rng_alloc(gsl_rng_rand48);
  unsigned int seed = PHRandomSeed();
  gsl_rng_set(fRandom,seed);
#endif
}
//=====
TPCSimulation::~TPCSimulation() {
  delete fPadMap;
  delete fEPDF;
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
  int mode=0;
  switch(mode) {
  case(0): // using analytical electron generation
    float ene = fEPDF->GetDEnergy();
    if(fVerbosity>20) std::cout << " Starting with " << ene << " keV" << std::endl;
#ifndef __CINT__
    ene = gsl_ran_poisson( fRandom, ene*kElectronsPerKeV );
#else
    ene = gRandom->PoissonD( ene*kElectronsPerKeV );
#endif
    if(fVerbosity>20) std::cout << " became " << ene << " electrons" << std::endl;
    while(ene>0) {
      float w = ene>kMaxElectronsPerHit?kMaxElectronsPerHit:ene;
      fEPDF->AddElectron(w);
      ene -= kMaxElectronsPerHit;
    }
    if(fVerbosity>20) std::cout << " Injected " << fEPDF->GetN() << " pdfs" << std::endl;

    break;
  }
  if(fVerbosity>2) std::cout << "Edep2Ele ... end." << std::endl;
}
//=====
void TPCSimulation::Transport()
{
  // compute path distortion due
  // to fields and ion backflow
  if(fVerbosity>2) std::cout << "Transport ... begin." << std::endl;
  int mode=0;
  switch(mode) {
  case(0): // using analytical symmetrical difussion and table for space charge
    float dlength = kGasHalfOfLength - TMath::Abs(fEPDF->GetZ());
    float tt = std::sqrt(kGasDiffusionTransverseSquare*dlength);
    float ll = std::sqrt(kGasDiffusionLongitudinalSquare*dlength);
    for(unsigned int i=0; i!=fEPDF->GetN(); ++i) {
#ifndef __CINT__
      float tra = gsl_ran_gaussian(fRandom,tt);
      float lon = gsl_ran_gaussian(fRandom,ll);
      float rph = gsl_ran_flat(fRandom,0,kTwoPi);
#else
      float tra = gRandom->Gaus(tt);
      float lon = gRandom->Gaus(ll);
      float rph = kTwoPi*gRandom->Rndm(1);
#endif
      float dx = tra*std::cos(rph);
      float dy = tra*std::sin(rph);
      float dz = lon;
      dx += 0.0;
      dy += 0.0;
      dz += 0.0;
      fEPDF->AddDx( i, dx );  // dx distortion
      fEPDF->AddDy( i, dy );  // dy distortion
      fEPDF->AddDz( i, dz );  // dz distortion
    }
    break;
  }
  if(fVerbosity>2) std::cout << "Transport ... end." << std::endl;
}
//====
void TPCSimulation::Amplify() {
  int mode=0;
  switch(mode) {
  case(0): // quick amplification
    for(unsigned int i=0; i!=fEPDF->GetN(); ++i) {
      fEPDF->Amplify( i, kAmplificationNominal );
      fEPDF->AddMST( i, kAmplificationSmearing );
    }
    break;
  }
}
//====
void TPCSimulation::Digitize() {
  int mode = 0;
  switch(mode) {
  case(0): // using two gaussian distributions: one for raise one for tail
    for(unsigned int i=0; i!=fEPDF->GetN(); ++i) {
      float dlength = kGasHalfOfLength-TMath::Abs(fEPDF->GetZ(i));
      fEPDF->SetT( i, dlength*kInverseOfDriftVelocity ); // map z(cm) to time
      fEPDF->SetMSL0( i, kTimeShapeRiseSquare );
      fEPDF->SetMSL1( i, kTimeShapeTailSquare );
    }
    break;
  }
}
//====
void TPCSimulation::PushToDigits(int blk) {
  for(unsigned int j=0; j!=fEPDF->GetN(); ++j) {
    if(fVerbosity>3) std::cout << "epdf #" << j << " || "
			       << " rad " << fEPDF->GetR(j) << " " << " phi " << fEPDF->GetPhi(j)
			       << " z " << fEPDF->GetZ(j) << " weight " << fEPDF->GetElectrons(j) << std::endl;
    ModuleRange_t mods = fPadMap->FindModuleRP( fEPDF->GetR(j), fEPDF->GetPhi(j),
						fEPDF->GetRMST(j), fEPDF->GetElectrons(j) );
    if(fVerbosity>2) std::cout << "range of modules " << mods.size() << std::endl;
    for(unsigned int k=0; k!=mods.size(); ++k)
      PushEPDF2Module(j,mods[k]+blk*kNModulesPerPlate);
  }
}
//====
void TPCSimulation::Hits2Digits(int blk) {
  unsigned int nhits = fHits->GetNHits(blk);
  if(fVerbosity>0) std::cout << "TPCSimulation::Hits2Digits(" << blk << ") => hits: " << nhits << std::endl;
  TPCHit *hit;
  for(unsigned int i=0; i!=nhits; ++i) {
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
    fEPDF->CopyFrom( hit );

    //======
    Edep2Ele();
    //======
    Transport();
    //======
    Amplify();
    //======
    Digitize();
    //======
    PushToDigits(blk);
  }
  if(fVerbosity>0) {
    int sum = 0;
    for(Module_t i=36*blk; i!=36+36*blk; ++i) sum += fDigits->GetNDigits(i);
    std::cout << "TPCSimulation::Hits2Digits(" << blk << ") => digits: " << sum << std::endl;
  }
}
//=====
void TPCSimulation::PushEPDF2Module(unsigned int clo, Module_t mod)
{
  if(fVerbosity>1) std::cout << "PushEPDF2Module " << int(mod) << std::endl;
  float rmst = fEPDF->GetRMST(clo);
  float cloR = fEPDF->GetR(clo);
  float cloP = fEPDF->GetPhi(clo);
  float cloT = fEPDF->GetT(clo);
  float rmsl0 = fEPDF->GetRMSL0(clo);
  float rmsl1 = fEPDF->GetRMSL1(clo);
  float nelec = fEPDF->GetElectrons(clo);
  if(fVerbosity>3) {
    std::cout << "CLOUD R " << cloR << " PHI " << cloP << " T " << cloT << " || "
	      << " RMST " << rmst << " RMSL " << rmsl0 << " " << rmsl1 << " || "
	      << " WEIGHT " << nelec << std::endl;
  }
  PadQuotaRange_t qRange = fPadMap->GetPadQuotasRP( nelec, mod, cloR, cloP, rmst );
  if(qRange.size()<1) return; // no coincidence found
  Pad_t prc;
  BinTime_t tim;
  float iww;
  Adc_t jww;
  if(fVerbosity>1) std::cout << "CLOUD WEIGHT " << nelec <<
		     " to be stored in " << qRange.size() <<  " pads." << std::endl;
  float sumiww = 0;
  for(unsigned int i=0; i!=qRange.size(); ++i) { // loop over traverse range
    prc = qRange[i].first;
    iww = qRange[i].second;
    sumiww += iww;
    TimeQuotaRange_t tRange = fPadMap->GetTimeQuotas( iww, cloT, rmsl0, rmsl1 );
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
