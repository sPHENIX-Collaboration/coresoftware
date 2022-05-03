// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_READEICFILES_H
#define G4MAIN_READEICFILES_H

#include <fun4all/SubsysReco.h>

#include <phhepmc/PHHepMCGenHelper.h>

#include <string>

class PHCompositeNode;
class TChain;

namespace erhic
{
class EventMC;
}

class ReadEICFiles : public SubsysReco, public PHHepMCGenHelper
{
 public:
  ReadEICFiles(const std::string &name = "EICReader");
  ~ReadEICFiles() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  /** Specify name of input file to open */
  bool OpenInputFile(const std::string &name);

  /** Set first entry from input tree to be used */
  void SetFirstEntry(int e) { entry = e; }
  /** Set name of output node */
  void SetNodeName(const std::string &s) { _node_name = s; }

 private:
  /** Get tree from input file */
  void GetTree();

  /** Creade node on node tree */
  int CreateNodeTree(PHCompositeNode *topNode);

  enum EvtGen
  {
    Milou = 1,
    DEMP = 2,
    Unknown = 100
  };

  /** Name of file containing input tree */
  std::string filename;

  /** Input tree created with eic-smear tree builder */
  TChain *Tin;

  /** Number of events in input tree */
  int nEntries;

  /** Number of current event being used from input tree */
  int entry;

  /** Event Generator id */
  int m_EvtGenId;

  /** Pinter to event record in tree (= branch).
      Use 'abstract' EventMC class pointer from which all
      event types (erhic::EventMilou etc) inherit from. */
  erhic::EventMC *GenEvent;

  // output
  std::string _node_name;

};

#endif /* READEICFILES_H__ */
