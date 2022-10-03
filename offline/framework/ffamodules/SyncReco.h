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

  void SegmentNumber(int i) {forced_segment = i;}

 private:
  int CreateNodeTree(PHCompositeNode *topNode);
// just if we need to override the segment for e.g. embedding
// where we want to reuse hijing files which normally set
// the segment number
  int forced_segment = -1;
};

#endif /* FFAMODULES_SYNCRECO_H */
