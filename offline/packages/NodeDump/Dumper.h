#ifndef NODEDUMP_DUMPER_H
#define NODEDUMP_DUMPER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHNodeDump;

class Dumper : public SubsysReco
{
 public:
  Dumper(const std::string &name = "DUMPER");
  virtual ~Dumper();
  int End(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void SetOutDir(const std::string &outdir);
  void SetPrecision(const int digits);
  int AddIgnore(const std::string &name);
  int Select(const std::string &name);

 private:
  PHNodeDump *nodedump;
};

#endif
