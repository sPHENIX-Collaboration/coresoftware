#ifndef BCOLUMICOUNT_STREAMINGBCOLUMICHECK_H
#define BCOLUMICOUNT_STREAMINGBCOLUMICHECK_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>

#include <TH1.h>


class StreamingBcoLumiCheck : public SubsysReco
{
 public:
  StreamingBcoLumiCheck(const std::string &name = "BCOLUMICHECKSTREAMINGOUTPUT");
  ~StreamingBcoLumiCheck() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
};

#endif  // BCOLUMICOUNT_STREAMINGBCOLUMICHECK_H
