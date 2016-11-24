#ifndef READEICFILES_H__
#define READEICFILES_H__

#include <fun4all/SubsysReco.h>

#include <string>

class PHComposteNode;
class TChain;

namespace erhic
{
  class EventMC;
}

class ReadEICFiles: public SubsysReco
{

 public:

  ReadEICFiles(const std::string &name="EICReader");
  virtual ~ReadEICFiles();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  bool OpenInputFile(const std::string &name);
  void GetTree();

  void SetFirstEntry(int e) {entry = e;}

 protected:

  TChain *Tin;
  /* use 'abstract' EventMC class pointer from which all
     event types (erhic::EventMilou etc) inherit from */
  erhic::EventMC * GenEvent;

  std::string filename;
  ///  Input Tree Variables
  int nEntries;
  int entry;
//  int ProcessID;
//  float Y;
//  float Q2;
//  float X;
//  float W2;
//  float NU;
};

#endif /* READEICFILES_H__ */
