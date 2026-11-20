#ifndef QA_TPCCALIBRATIONQA_RUNAVERAGE_H
#define QA_TPCCALIBRATIONQA_RUNAVERAGE_H
#include <fun4all/SubsysReco.h>
class PHCompositeNode;
class TpcRunAverageConditions : public SubsysReco
{
 public:
  TpcRunAverageConditions();
  int InitRun(PHCompositeNode *topNode) override;
};
#endif
