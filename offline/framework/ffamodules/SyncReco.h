#ifndef FFAMODULES_SYNCRECO_H
#define FFAMODULES_SYNCRECO_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SyncReco : public SubsysReco
{
 public:
  SyncReco(const std::string &name = "SYNC");
  ~SyncReco() override {}

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 protected:
  int CreateNodeTree(PHCompositeNode *topNode);
};

#endif /* FFAMODULES_SYNCRECO_H */
