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

  ~SpinDBNode() override = default;
  int InitRun(PHCompositeNode *topNode) override;
  void SetQALevel(int qa_level) { m_qa_level = qa_level; }

 private:
    int m_qa_level{-999};
};

#endif  // SPINDBNODE_H
