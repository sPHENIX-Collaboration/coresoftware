#ifndef TRACKBASE_TRKRTRUTHTRACKV1_H
#define TRACKBASE_TRKRTRUTHTRACKV1_H
/**
 * @file g4tpc/TrkrTruthTrackv1.h
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of TrkrTruthTrack
 */
#include "TrkrTruthTrack.h"

class PHG4Particle;
class PHG4VtxPoint;

class TrkrTruthTrackv1 : public TrkrTruthTrack
{
 public:
  //! ctor
  TrkrTruthTrackv1();
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

  void setTrackid(unsigned int _) override { trackid = _; };

  void setX0(float _) override { X0 = _; };
  void setY0(float _) override { Y0 = _; };
  void setZ0(float _) override { Z0 = _; };

  void setPseudoRapity(float _) override { pseudoRapidity = _; };
  void setPt(float _) override { pt = _; };
  void setPhi(float _) override { phi = _; };

  void identify(std::ostream& os = std::cout) const override;

 protected:
  unsigned int trackid;
  float X0;
  float Y0;
  float Z0;

  float pseudoRapidity;
  float pt;
  float phi;

  std::vector<TrkrDefs::cluskey> clusters;
  ClassDefOverride(TrkrTruthTrackv1, 1)
};
#endif
