#ifndef HEPMCCOMPRESS_H__
#define HEPMCCOMPRESS_H__

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;

class HepMCCompress : public SubsysReco
{
 public:
  HepMCCompress(const std::string &name = "HEPMCREADER");
  virtual ~HepMCCompress() {}

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
 protected:
  short int FloatToInt(const float rval) const;
  std::set<int> exclude_pid;
  std::set<int> select_pid;
};

#endif

