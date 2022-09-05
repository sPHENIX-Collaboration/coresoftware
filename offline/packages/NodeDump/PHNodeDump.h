#ifndef NODEDUMP_PHNODEDUMP_H
#define NODEDUMP_PHNODEDUMP_H

#include <phool/PHNodeOperation.h>

#include <map>
#include <set>
#include <string>

class PHNode;
class DumpObject;
class PHCompositeNode;

class PHNodeDump : public PHNodeOperation
{
 public:
  PHNodeDump() {}
  ~PHNodeDump() override;
  int CloseOutputFiles();
  int GetGlobalVars(PHCompositeNode *topNode);

  int RunNumber() { return runnumber; }
  int EvtSequence() { return evtsequence; }

  int AddIgnore(const std::string &name);
  int Select(const std::string &name);
  void SetOutDir(const std::string &dirname) { outdir = dirname; }
  void SetPrecision(const int digits) { fp_precision = digits; }
  void PrintEvtSeq(const int i) { print_evtseq = i; }

 private:
  void perform(PHNode *) override;
  int AddDumpObject(const std::string &NodeName, PHNode *node);
  std::map<std::string, DumpObject *> dumpthis;
  std::set<std::string> ignore;
  std::set<std::string> exclusive;
  int runnumber = 0;
  int evtsequence = -9999;
  int fp_precision = -1;
  std::string outdir = "./";
  int print_evtseq = 1;
};

#endif
