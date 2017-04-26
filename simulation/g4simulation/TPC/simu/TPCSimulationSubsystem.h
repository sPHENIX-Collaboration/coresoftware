#ifndef __TPCSIMULATIONSUBSYSTEM_H__
#define __TPCSIMULATIONSUBSYSTEM_H__

#include <fun4all/SubsysReco.h>
#include <TString.h>

class PHCompositeNode;
class TPCSimulation;
class TStopwatch;
class TH1F;

class TPCSimulationSubsystem : public SubsysReco {
 public:
  TPCSimulationSubsystem();
  virtual ~TPCSimulationSubsystem();
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void SetTreeFileName(TString name) {fTreeFileName=name;}
  
 protected:
  TPCSimulation *fSimulation;
  TStopwatch *fStopwatch;
  TH1F *fHist;

  TString fTreeFileName;
};

#endif
