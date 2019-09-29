#ifndef __TPCIntegratedCharge_H__
#define __TPCIntegratedCharge_H__

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Fun4AllHistoManager;

/// \class TPCIntegratedCharge
class TPCIntegratedCharge : public SubsysReco
{
 public:
  TPCIntegratedCharge(
      unsigned int minLayer,
      unsigned int m_maxLayer,
      const std::string &outputfilename = "TPCIntegratedCharge.root");

  virtual ~TPCIntegratedCharge();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
#if !defined(__CINT__) || defined(__CLING__)

  Fun4AllHistoManager *getHistoManager();

  std::string m_outputFileName;

  unsigned int m_minLayer;
  unsigned int m_maxLayer;

#endif  // #if !defined(__CINT__) || defined(__CLING__)
};

#endif  // __CALOEVALUATOR_H__
