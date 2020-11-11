
#include <Acts/Seeding/InternalSpacePoint.hpp>
#include <Acts/Seeding/InternalSeed.hpp>
#include <Acts/Seeding/IExperimentCuts.hpp>

#include <algorithm>

template <typename SpacePoint>
class SphenixActsDetectorCuts : public Acts::IExperimentCuts<SpacePoint> {
 public:
  struct Config {
    int maxSeedSize = 5;
    float cutRadius = 200.;
    float cutWeight = 380.;
    float weight_outer = 400.;
    float weight_inner = 200.;
    float keepWeight = 200.;
    float minRadius = 2.;
  };

  /// Returns seed weight bonus/malus depending on detector considerations.
  float seedWeight(const Acts::InternalSpacePoint<SpacePoint>& bottom,
                   const Acts::InternalSpacePoint<SpacePoint>& middle,
                   const Acts::InternalSpacePoint<SpacePoint>& top) const;
 
  bool singleSeedCut(float weight, const Acts::InternalSpacePoint<SpacePoint>& bottom,
                     const Acts::InternalSpacePoint<SpacePoint>& middle,
                     const Acts::InternalSpacePoint<SpacePoint>& top) const;

  std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<SpacePoint>>>>
  cutPerMiddleSP(
      std::vector<
          std::pair<float, std::unique_ptr<const Acts::InternalSeed<SpacePoint>>>>
          seeds) const;

 private:
  Config m_cfg;
};

template <typename SpacePoint>
float SphenixActsDetectorCuts<SpacePoint>::seedWeight(
    const Acts::InternalSpacePoint<SpacePoint>& bottom,
    const Acts::InternalSpacePoint<SpacePoint>&,
    const Acts::InternalSpacePoint<SpacePoint>& top) const {
  float weight = 0;
  if (bottom.radius() > m_cfg.cutRadius) {
    weight = m_cfg.weight_outer;
  }
  if (top.radius() < m_cfg.cutRadius) {
    weight = m_cfg.weight_inner;
  }
  return weight;
}

template <typename SpacePoint>
bool SphenixActsDetectorCuts<SpacePoint>::singleSeedCut(
    float weight, const Acts::InternalSpacePoint<SpacePoint>& b,
    const Acts::InternalSpacePoint<SpacePoint>&,
    const Acts::InternalSpacePoint<SpacePoint>&) const {
  return !(b.radius() > m_cfg.cutRadius && weight < m_cfg.cutWeight);
}

template <typename SpacePoint>
std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<SpacePoint>>>>
SphenixActsDetectorCuts<SpacePoint>::cutPerMiddleSP(
    std::vector<
        std::pair<float, std::unique_ptr<const Acts::InternalSeed<SpacePoint>>>>
        seeds) const {
  std::vector<std::pair<float, std::unique_ptr<const Acts::InternalSeed<SpacePoint>>>>
      newSeedsVector;
  if (seeds.size() > 1) {
    newSeedsVector.push_back(std::move(seeds[0]));
    size_t itLength = std::min(seeds.size(), size_t(m_cfg.maxSeedSize));
    // don't cut first element
    for (size_t i = 1; i < itLength; i++) {
      if (seeds[i].first > m_cfg.keepWeight ||
          seeds[i].second->sp[0]->radius() > m_cfg.minRadius) {
        newSeedsVector.push_back(std::move(seeds[i]));
      }
    }
    return newSeedsVector;
  }
  return seeds;
}
