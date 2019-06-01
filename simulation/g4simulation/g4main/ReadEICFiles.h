// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_READEICFILES_H
#define G4MAIN_READEICFILES_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <string>

class PHHepMCGenEvent;
class TChain;

namespace erhic
{
class EventMC;
}

class ReadEICFiles : public SubsysReco
{
 public:
  ReadEICFiles(const std::string &name = "EICReader");
  virtual ~ReadEICFiles();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  /** Specify name of input file to open */
  bool OpenInputFile(const std::string &name);

  /** Set first entry from input tree to be used */
  void SetFirstEntry(int e) { entry = e; }
  /** Set name of output node */
  void SetNodeName(const std::string &s) { _node_name = s; }
  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(PHHepMCGenHelper::VTXFUNC x, PHHepMCGenHelper::VTXFUNC y, PHHepMCGenHelper::VTXFUNC z, PHHepMCGenHelper::VTXFUNC t)
  {
    hepmc_helper.set_vertex_distribution_function(x, y, z, t);
  }

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_mean(x, y, z, t);
  }

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_width(x, y, z, t);
  }
  //
  //! reuse vertex from another PHHepMCGenEvent with embedding_id = src_embedding_id Additional smearing and shift possible with set_vertex_distribution_*()
  void set_reuse_vertex(int src_embedding_id)
  {
    hepmc_helper.set_reuse_vertex(src_embedding_id);
  }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return hepmc_helper.get_embedding_id(); }
  //
  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { hepmc_helper.set_embedding_id(id); }
 protected:
  /** Get tree from input file */
  void GetTree();

  /** Creade node on node tree */
  int CreateNodeTree(PHCompositeNode *topNode);

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
  erhic::EventMC *GenEvent;

  // output
  std::string _node_name;

  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper hepmc_helper;
};

#endif /* READEICFILES_H__ */
