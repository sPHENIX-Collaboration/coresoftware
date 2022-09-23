// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SVTXTRACKSTATEREMOVAL_H
#define SVTXTRACKSTATEREMOVAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SvtxTrackStateRemoval : public SubsysReco
{
 public:
  SvtxTrackStateRemoval(const std::string &name = "SvtxTrackStateRemoval");

  ~SvtxTrackStateRemoval() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
};

#endif  // SVTXTRACKSTATEREMOVAL_H
