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

  /** Specify name of input file to open */
  bool OpenInputFile(const std::string &name);

  /** Set first entry from input tree to be used */
  void SetFirstEntry(int e) {entry = e;}

 protected:

  /** Get tree from input file */
  void GetTree();

  /** Name of file containing input tree */
  std::string filename;

  /** Input tree created with eic-smear tree builder */
  TChain *Tin;

  /** Number of events in input tree */
  int nEntries;

  /** Number of current event being used from input tree */
  int entry;

  /** Pinter to event record in tree (= branch).
      Use 'abstract' EventMC class pointer from which all
      event types (erhic::EventMilou etc) inherit from. */
  erhic::EventMC * GenEvent;

};

#endif /* READEICFILES_H__ */
