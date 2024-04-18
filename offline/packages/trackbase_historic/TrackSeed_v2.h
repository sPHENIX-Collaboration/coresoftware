#ifndef TRACKBASEHISTORIC_TRACKSEED_V2_H
#define TRACKBASEHISTORIC_TRACKSEED_V2_H

#include "TrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <limits.h>
#include <cmath>
#include <iostream>

class TrackSeed_v2 : public TrackSeed
{
 public:
  TrackSeed_v2();

  /// Copy constructors
  TrackSeed_v2(const TrackSeed&);
  TrackSeed_v2(const TrackSeed_v2&);
  TrackSeed_v2& operator=(const TrackSeed_v2& seed);
  ~TrackSeed_v2() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrackSeed_v2(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed&) override;
  void CopyFrom(TrackSeed* seed) override { CopyFrom(*seed); }
  PHObject* CloneMe() const override { return new TrackSeed_v2(*this); }

  // method to return phi from a given set of global positions
  float get_phi(const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions) const override;   // returns phi calculated from supplied cluster positions

  // methods that return values based on track fit parameters
  float get_pz() const override;
  float get_x() const override;
  float get_y() const override;
  float get_z() const override;
  float get_eta() const override;
  float get_theta() const override;
  float get_pt() const override;
  float get_p() const override;
  float get_px() const override;
  float get_py() const override;

  //methods that return member variables
  int get_charge() const override;
  float get_qOverR() const override { return m_qOverR; }
  float get_X0() const override { return m_X0; }
  float get_Y0() const override { return m_Y0; }
  float get_slope() const override { return m_slope; }
  float get_Z0() const override { return m_Z0; }
  float get_phi() const override  { return m_phi; }  // returns the stored phi
  short int get_crossing() const override { return m_crossing; }  
  void set_crossing(const short int crossing) override { m_crossing = crossing; }
  void set_qOverR(const float qOverR) override { m_qOverR = qOverR; }
  void set_X0(const float X0) override { m_X0 = X0; }
  void set_Y0(const float Y0) override { m_Y0 = Y0; }
  void set_slope(const float slope) override { m_slope = slope; }
  void set_Z0(const float Z0) override { m_Z0 = Z0; }
  void set_phi(const float phi) override { m_phi = phi; }

  void clear_cluster_keys() override { m_cluster_keys.clear(); }
  bool empty_cluster_keys() const override { return m_cluster_keys.empty(); }
  size_t size_cluster_keys() const override { return m_cluster_keys.size(); }

  void insert_cluster_key(TrkrDefs::cluskey clusterid) override { m_cluster_keys.insert(clusterid); }
  size_t erase_cluster_key(TrkrDefs::cluskey clusterid) override { return m_cluster_keys.erase(clusterid); }
  ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const override { return m_cluster_keys.find(clusterid); }
  ConstClusterKeyIter begin_cluster_keys() const override { return m_cluster_keys.begin(); }
  ConstClusterKeyIter end_cluster_keys() const override { return m_cluster_keys.end(); }
  ClusterKeyIter find_cluster_keys(unsigned int clusterid) override { return m_cluster_keys.find(clusterid); }
  ClusterKeyIter begin_cluster_keys() override { return m_cluster_keys.begin(); }
  ClusterKeyIter end_cluster_keys() override { return m_cluster_keys.end(); }

  /// Updates R, X0, Y0
  void circleFitByTaubin(const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
                         uint8_t startLayer = 0,
                         uint8_t endLayer = 58) override;

  /// Updates r-z slope and intercept B
  void lineFit(const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
               uint8_t startLayer = 0,
               uint8_t endLayer = 58) override;

 protected:
  /// Returns transverse PCA to (0,0)
  std::pair<float, float> findRoot() const;

 private:
  ClusterKeySet m_cluster_keys;

  float m_qOverR = NAN;
  float m_X0 = NAN;
  float m_Y0 = NAN;
  float m_slope = NAN;
  float m_Z0 = NAN;
  float m_phi = NAN;

  short int m_crossing = std::numeric_limits<short int>::max();

  ClassDefOverride(TrackSeed_v2, 1);
};

#endif
