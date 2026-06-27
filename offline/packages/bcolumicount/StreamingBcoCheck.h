#ifndef BCOLUMICOUNT_STREAMINGBCOCHECK_H
#define BCOLUMICOUNT_STREAMINGBCOCHECK_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <string>

#include <TH1.h>


class StreamingBcoCheck : public SubsysReco
{
 public:
  StreamingBcoCheck(const std::string &name = "BCOCHECKSTREAMINGOUTPUT");
  ~StreamingBcoCheck() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
};

#endif  // BCOLUMICOUNT_STREAMINGBCOCHECK_H
