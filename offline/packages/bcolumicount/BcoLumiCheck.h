#ifndef BCOLUMICOUNT_BCOLUMICHECK_H
#define BCOLUMICOUNT_BCOLUMICHECK_H

#include <fun4all/SubsysReco.h>

#include <string>

class BcoLumiCheck : public SubsysReco
{
 public:
  BcoLumiCheck(const std::string &name = "BCOLUMICHECK");
  ~BcoLumiCheck() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
};

#endif  // BCOLUMICOUNT_BCOLUMICHECK_H
