// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef SPINDBNODE_H
#define SPINDBNODE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SpinDBNode : public SubsysReco
{
 public:

  SpinDBNode(const std::string &name = "SpinDBNode");

  ~SpinDBNode() override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
};

#endif // SPINDBNODE_H
