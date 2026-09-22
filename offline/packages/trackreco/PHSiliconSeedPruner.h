// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef PHSILICONSEEDPRUNER_H
#define PHSILICONSEEDPRUNER_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <cstddef>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;

class PHSiliconSeedPruner : public SubsysReco
{
 public:
  PHSiliconSeedPruner(const std::string& name = "PHSiliconSeedPruner");
  ~PHSiliconSeedPruner() override;

  PHSiliconSeedPruner(const PHSiliconSeedPruner&) = delete;
  PHSiliconSeedPruner& operator=(const PHSiliconSeedPruner&) = delete;
  PHSiliconSeedPruner(PHSiliconSeedPruner&&) = delete;
  PHSiliconSeedPruner& operator=(PHSiliconSeedPruner&&) = delete;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void set_track_map_name(const std::string& name) { m_trackMapName = name; }
  void set_output_track_map_name(const std::string& name) { m_outputTrackMapName = name; }
  void set_random_seed(unsigned int seed);

 private:
  static constexpr std::size_t kMvtxLayerCount = 3;
  static constexpr std::size_t kMaxRepresentatives = 10;
  static constexpr std::size_t kSearchBudget = 2000000;
  static constexpr std::size_t kHeuristicRestarts = 256;
  static constexpr std::size_t kDebugMinimumGroupSize = 10;

  TrackSeedContainer* m_siliconSeeds{nullptr};
  TrackSeedContainer* m_prunedSiliconSeeds{nullptr};
  std::string m_trackMapName{"SiliconTrackSeedContainer"};              // the default input node name for the TrackSeedContainer
  std::string m_outputTrackMapName{"PrunedSiliconTrackSeedContainer"};  // the default output node name for the TrackSeedContainer

  unsigned int m_randomSeed{0};  // the default random seed
  gsl_rng* m_rng{nullptr};
};

#endif  // PHSILICONSEEDPRUNER_H
