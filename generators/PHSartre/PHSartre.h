#ifndef PHSARTRE_PHSARTRE_H
#define PHSARTRE_PHSARTRE_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <cmath>
#include <string>
#include <vector>

class PHCompositeNode;

class Sartre;
class Event;
class EventGeneratorSettings;
class PHSartreGenTrigger;
class TGenPhaseSpace;

class PHSartre : public SubsysReco
{
 public:
  PHSartre(const std::string &name = "PHSartre");
  virtual ~PHSartre();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_config_file(const std::string &cfg_file)
  {
    _configFile = cfg_file;
  }

  void print_config() const;

  /// set event selection criteria
  void register_trigger(PHSartreGenTrigger *theTrigger);
  void set_trigger_OR()
  {
    _triggersOR = true;
    _triggersAND = false;
  }  // default true
  void set_trigger_AND()
  {
    _triggersAND = true;
    _triggersOR = false;
  }

  /// pass commands directly to PYTHIA8
  void process_string(std::string s) { _commands.push_back(s); }

  void beam_vertex_parameters(double beamX,
                              double beamY,
                              double beamZ,
                              double beamXsigma,
                              double beamYsigma,
                              double beamZsigma)
  {
    set_vertex_distribution_mean(beamX, beamY, beamZ, 0);
    set_vertex_distribution_width(beamXsigma, beamYsigma, beamZsigma, 0);
  }

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

 private:
  int create_node_tree(PHCompositeNode *topNode);
  double percent_diff(const double a, const double b) { return fabs((a - b) / a); }
  void randomlyReverseBeams(Event *myEvent);
  void ReverseBeams(Event *myEvent);

  int _eventcount;
  int _gencount;

  // event selection
  std::vector<PHSartreGenTrigger *> _registeredTriggers;
  bool _triggersOR;
  bool _triggersAND;

  std::string _configFile;
  std::vector<std::string> _commands;

  // Sartre
  Sartre *_sartre;
  EventGeneratorSettings *settings;
  TGenPhaseSpace *decay;
  int daughterID;
  double daughterMasses[2];
  bool doPerformDecay;

  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper hepmc_helper;
};

#endif /* __PHSARTRE_H__ */
