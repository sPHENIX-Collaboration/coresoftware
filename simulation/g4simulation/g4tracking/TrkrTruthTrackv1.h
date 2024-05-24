#ifndef TRACKBASE_TRKRTRUTHTRACKV1_H
#define TRACKBASE_TRKRTRUTHTRACKV1_H
/**
 * @file g4tpc/TrkrTruthTrackv1.h
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of TrkrTruthTrack
 */
#include "TrkrTruthTrack.h"

#include <limits>

class PHG4Particle;
class PHG4VtxPoint;

class TrkrTruthTrackv1 : public TrkrTruthTrack
{
 public:
  //! ctor
  TrkrTruthTrackv1() = default;
  TrkrTruthTrackv1(unsigned int, PHG4Particle*, PHG4VtxPoint*);

  //! dtor
  ~TrkrTruthTrackv1() override = default;

  unsigned int getTrackid() const override { return trackid; };
  float getX0() const override { return X0; };
  float getY0() const override { return Y0; };
  float getZ0() const override { return Z0; };

  float getPseudoRapidity() const override { return pseudoRapidity; };
  float getPt() const override { return pt; };
  float getPhi() const override { return phi; };

  std::vector<TrkrDefs::cluskey>& getClusters() override { return clusters; };
  void addCluster(TrkrDefs::cluskey) override;

  bool has_hitsetkey(TrkrDefs::hitsetkey) const override;  // note, only works when the data is sorted
  bool has_hitsetkey(TrkrDefs::cluskey) const override;
  std::pair<bool, TrkrDefs::cluskey> get_cluskey(TrkrDefs::hitsetkey) const override;  // bool is if there is the key, and

  void setTrackid(unsigned int _) override { trackid = _; };

  void setX0(float _) override { X0 = _; };
  void setY0(float _) override { Y0 = _; };
  void setZ0(float _) override { Z0 = _; };

  void setPseudoRapity(float _) override { pseudoRapidity = _; };
  void setPt(float _) override { pt = _; };
  void setPhi(float _) override { phi = _; };

  void identify(std::ostream& os = std::cout) const override;

  struct CompHitSetKey
  {
    bool operator()(const TrkrDefs::cluskey& lhs, const TrkrDefs::cluskey& rhs)
    {
      return TrkrDefs::getHitSetKeyFromClusKey(lhs) < TrkrDefs::getHitSetKeyFromClusKey(rhs);
    };
    bool operator()(const TrkrDefs::cluskey& lhs, const TrkrDefs::hitsetkey& rhs)
    {
      return TrkrDefs::getHitSetKeyFromClusKey(lhs) < rhs;
    };
    bool operator()(const TrkrDefs::hitsetkey& lhs, const TrkrDefs::cluskey& rhs)
    {
      return lhs < TrkrDefs::getHitSetKeyFromClusKey(rhs);
    };
  };

 protected:
  unsigned int trackid{std::numeric_limits<unsigned int>::max()};
  float X0{std::numeric_limits<float>::quiet_NaN()};
  float Y0{std::numeric_limits<float>::quiet_NaN()};
  float Z0{std::numeric_limits<float>::quiet_NaN()};

  float pseudoRapidity{std::numeric_limits<float>::quiet_NaN()};
  float pt{std::numeric_limits<float>::quiet_NaN()};
  float phi{std::numeric_limits<float>::quiet_NaN()};

  std::vector<TrkrDefs::cluskey> clusters;
  ClassDefOverride(TrkrTruthTrackv1, 1)
};
#endif
