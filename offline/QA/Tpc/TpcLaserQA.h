#ifndef QA_TPC_TPCLASERQA_H
#define QA_TPC_TPCLASERQA_H


#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class LaserClusterContainer;

class TH1;
class TH2;

class TpcLaserQA : public SubsysReco
{
	
public:
  explicit TpcLaserQA(const std::string &name = "TpcLaserQA");
  
  ~TpcLaserQA() override = default;
  
  int InitRun(PHCompositeNode *topNode) override;  
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  
private:

  void createHistos();
  std::string getHistoPrefix() const;
  
  double rBinEdges[5] = {0.0, 0.256, 0.504, 0.752, 1.0};

  TH1* m_nLaserEvents{nullptr};
  TH2* m_TPCWheel[2]{nullptr};

};

#endif
