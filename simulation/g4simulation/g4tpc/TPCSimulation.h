#ifndef __TPCSIMULATION_H__
#define __TPCSIMULATION_H__

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <TPCDataTypes.h>
#include <TString.h>
using namespace TPCDataTypes;

class TFile;
class TTree;
class TPCElectronPDF;
class TPCHitsContainer;
class TPCDigitsContainer;
class TPCPadMap;

class TPCSimulation {
 public:
  TPCSimulation();
  virtual ~TPCSimulation();
  void ConnectHits(TPCHitsContainer *in) {fHits = in;}
  void ConnectDigits(TPCDigitsContainer *out) {fDigits = out;}
  void Hits2Digits();
  void Hits2Digits(int);
  void SetVerbosity(int v) {fVerbosity=v;}
  void PrepareTree(TString name);
  void WriteFile();
  
 protected:
  void Edep2Ele();
  void Transport();
  void Amplify();
  void Digitize();
  void PushToDigits(int);
  void PushEPDF2Module(unsigned int i, Module_t mod);
  TimeQuotaRange_t MakeTQuotas( Float_t ele, Float_t dlength, Float_t rms);
  
#ifndef __CINT__
  gsl_rng *fRandom;
#endif

  TPCHitsContainer *fHits; ///>! container for tpchits
  TPCDigitsContainer *fDigits; ///>! container for tpcdigits
  TPCPadMap *fPadMap;
  TPCElectronPDF *fEPDF;
  int fVerbosity;
  TString fTreeFileName;
  TTree *fTree;
  TFile *fFile;
  typedef struct {
    Int_t track;
    Float_t r;
    Float_t phi;
    Float_t z;
    Float_t length;
    Float_t edep;
  } TreeHitRegister_t;
  TreeHitRegister_t fHitRegister;
};

#endif
