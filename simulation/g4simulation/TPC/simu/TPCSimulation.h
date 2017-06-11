#ifndef __TPCSIMULATION_H__
#define __TPCSIMULATION_H__

#include <TPCbase/TPCDataTypes.h>
#include <TString.h>
using namespace TPCDataTypes;

class TFile;
class TTree;
class TPCCloud;
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
  void PushCloud2Module(Module_t mod);
  TimeQuotaRange_t MakeTQuotas( Float_t ele, Float_t dlength, Float_t rms);

  TPCHitsContainer *fHits; ///>! container for tpchits
  TPCDigitsContainer *fDigits; ///>! container for tpcdigits
  TPCPadMap *fPadMap;
  TPCCloud *fCloud;
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
