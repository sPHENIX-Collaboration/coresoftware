#ifndef FFAMODULES_SYNCRECO_H
#define FFAMODULES_SYNCRECO_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SyncReco : public SubsysReco
{
 public:
  SyncReco(const std::string &name = "SYNC");
  virtual ~SyncReco() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 protected:
  int CreateNodeTree(PHCompositeNode *topNode);
};

#endif /* FFAMODULES_SYNCRECO_H */
