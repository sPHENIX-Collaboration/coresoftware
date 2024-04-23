// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_GL1CHECK_H
#define FFARAWMODULES_GL1CHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class Gl1Check : public SubsysReco, public DumpPacket
{
 public:
  Gl1Check(const std::string &name = "Gl1Check");

  ~Gl1Check() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  //  int ResetEvent(PHCompositeNode *topNode) override;

  //  void MyEvtNode(const std::string &name) {m_EvtNodeName = name;}

 private:
  /* std::string m_EvtNodeName = "GL1RAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_GL1CHECK_H
