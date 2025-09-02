#ifndef TRUTHNEUTRALMESONBUILDER_H
#define TRUTHNEUTRALMESONBUILDER_H

#include <fun4all/SubsysReco.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4VtxPoint.h>

#include <map>
#include <set>
#include <string>
#include <unordered_map>

class PHCompositeNode;
class PHG4Shower;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class TruthNeutralMesonContainer;
namespace HepMC
{
  class GenParticle;
}

class TruthNeutralMesonBuilder : public SubsysReco
{
 public:
  explicit TruthNeutralMesonBuilder(const std::string &name = "TruthNeutralMesonBuilder");
  ~TruthNeutralMesonBuilder() override = default;

  int Init(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void classify_eta_decay(int bc, HepMC::GenParticle *hepMom, bool &eta_3pi0, bool &eta_pi0pipm);

  static bool FindConversion(PHG4TruthInfoContainer *_truthinfo, int trackid, float energy);

  int process_event(PHCompositeNode *topNode) override;

  void set_truthnmeson_NodeName(const std::string &inputNodeName) { m_truthnmeson_node_name = inputNodeName; }

  void set_save_eta(bool _sv) { m_save_eta = _sv; }
  void set_save_pi0(bool _sv) { m_save_pi0 = _sv; }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  bool m_save_eta{true};
  bool m_save_pi0{true};

  std::string m_truthnmeson_node_name = "TruthNeutralMeson";

  std::set<std::pair<int, int>> m_seen_mother_keys;

  PHG4TruthInfoContainer *truthinfo{nullptr};
  PHHepMCGenEventMap *genevtmap{nullptr};
  TruthNeutralMesonContainer *_container{nullptr};

  std::unordered_map<int, PHG4Particle *> g4_by_barcode;
  std::unordered_map<int, PHG4Particle *> g4_by_id;
  std::unordered_map<int, std::vector<PHG4Particle *>> g4_children;
};

#endif
