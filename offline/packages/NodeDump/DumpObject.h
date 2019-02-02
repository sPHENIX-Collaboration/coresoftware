#ifndef NODEDUMP_DUMPOBJECT_H
#define NODEDUMP_DUMPOBJECT_H

#include <string>
#include <iosfwd>

class PHNode;
class PHNodeDump;

class DumpObject
{
 public:
  DumpObject(const std::string &NodeName = "DUMMY");
  virtual ~DumpObject() {}

  virtual int Init(); // called during intialization
  virtual int process_event(PHNode *mynode);
  virtual int CloseOutputFile();
  virtual const char *Name() const
    {
      return ThisName.c_str();
    }
  virtual void Print(const char *what) const;
  virtual void SetOutDir(const std::string &outdir) {OutDir = outdir;}
  virtual void SetPrecision(const int digits) {fp_precision = digits;}
  void SetParentNodeDump(PHNodeDump *nd) {myNodeDump = nd;}
  void NoOutput(const int i=1) {no_output = i;}
  void WriteRunEvent(int i) {write_run_event = i;}

 protected:
  virtual int process_Node(PHNode *myNode);
  virtual int OpenOutFile();
  std::ostream *fout;

 private:
  std::string ThisName;
  int verbosity;
  int write_run_event; // flag for not writing info for each event
  std::string OutDir;
  int fp_precision;
  PHNodeDump *myNodeDump;
  int no_output;
};

#endif

