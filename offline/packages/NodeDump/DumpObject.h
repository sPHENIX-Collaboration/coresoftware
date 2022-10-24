#ifndef NODEDUMP_DUMPOBJECT_H
#define NODEDUMP_DUMPOBJECT_H

#include <iosfwd>
#include <string>

class PHNode;
class PHNodeDump;

class DumpObject
{
 public:
  DumpObject(const std::string &NodeName = "DUMMY");
  virtual ~DumpObject() {}

  virtual int Init();  // called during intialization
  virtual int process_event(PHNode *mynode);
  virtual int CloseOutputFile();
  virtual const std::string Name() const { return ThisName; }
  virtual void Print(const char *what) const;
  virtual void SetOutDir(const std::string &outdir) { OutDir = outdir; }
  virtual void SetPrecision(const int digits) { fp_precision = digits; }
  void SetParentNodeDump(PHNodeDump *nd) { myNodeDump = nd; }
  void NoOutput(const int i = 1) { no_output = i; }
  void WriteRunEvent(const int i) { write_run_event = i; }
  void PrintEvtSeq(const int i) { print_evtseq = i; }

 protected:
  virtual int process_Node(PHNode *myNode);
  virtual int OpenOutFile();
  std::ostream *fout = nullptr;

 private:
  std::string ThisName;
  int write_run_event = 1;  // flag for not writing info for each event
  std::string OutDir = "./";
  int fp_precision = -1;
  PHNodeDump *myNodeDump = nullptr;
  int no_output = 0;
  int print_evtseq = 1;
};

#endif
