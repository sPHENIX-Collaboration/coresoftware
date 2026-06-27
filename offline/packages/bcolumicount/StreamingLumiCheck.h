#ifndef BCOLUMICOUNT_STREAMINGLUMICHECK_H
#define BCOLUMICOUNT_STREAMINGLUMICHECK_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>

#include <TH1.h>


class StreamingLumiCheck : public SubsysReco
{
 public:
  StreamingLumiCheck(const std::string &name = "LUMICHECKSTREAMINGOUTPUT");
  ~StreamingLumiCheck() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  //int process_event(PHCompositeNode *topNode) override;

 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
};

#endif  // BCOLUMICOUNT_STREAMINGLUMICHECK_H
