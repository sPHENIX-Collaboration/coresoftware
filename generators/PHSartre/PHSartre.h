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

class PHSartre : public SubsysReco, public PHHepMCGenHelper
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

 private:
  double percent_diff(const double a, const double b) { return fabs((a - b) / a); }
  void randomlyReverseBeams(Event *myEvent);
  void ReverseBeams(Event *myEvent);

  int _eventcount = 0;
  int _gencount = 0;

  // event selection
  std::vector<PHSartreGenTrigger *> _registeredTriggers;
  bool _triggersOR = true;
  bool _triggersAND = false;

  std::string _configFile;
  std::vector<std::string> _commands;

  // Sartre
  Sartre *_sartre = nullptr;
  EventGeneratorSettings *settings = nullptr;
  TGenPhaseSpace *decay = nullptr;
  int daughterID = -1;
  double daughterMasses[2] = {};
  bool doPerformDecay = false;
};

#endif /* __PHSARTRE_H__ */
