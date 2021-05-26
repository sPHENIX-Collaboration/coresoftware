// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_HEPMCCOMPRESS_H
#define G4MAIN_HEPMCCOMPRESS_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;

class HepMCCompress : public SubsysReco
{
 public:
  HepMCCompress(const std::string &name = "HEPMCREADER");
  ~HepMCCompress() override {}

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
 protected:
  short int FloatToInt(const float rval) const;
  std::set<int> exclude_pid;
  std::set<int> select_pid;
};

#endif

