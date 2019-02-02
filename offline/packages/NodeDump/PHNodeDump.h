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
  PHNodeDump(); 
  virtual ~PHNodeDump();
  int CloseOutputFiles();
  int GetGlobalVars(PHCompositeNode *topNode);

  int RunNumber() {return runnumber;}
  int EvtSequence() {return evtsequence;}

  int AddIgnore(const std::string &name);
  int Select(const std::string &name);
  int SetOutDir(const std::string &dirname);
  void SetPrecision(const int digits) {fp_precision = digits;}

protected: 
   virtual void perform(PHNode*);
   int AddDumpObject(const std::string &NodeName, PHNode *node);
   std::map <std::string, DumpObject *> dumpthis;
   std::set<std::string> ignore;
   std::set<std::string> exclusive;
   int runnumber;
   int evtsequence;
   int fp_precision;
   std::string outdir;
}; 

#endif
