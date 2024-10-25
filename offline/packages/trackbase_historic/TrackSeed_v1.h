#ifndef TRACKBASEHISTORIC_TRACKSEED_V1_H
#define TRACKBASEHISTORIC_TRACKSEED_V1_H

#include "TrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <limits.h>
#include <cmath>
#include <iostream>

class TrackSeed_v1 : public TrackSeed
{
 public:
  TrackSeed_v1() = default;

  /// Copy constructors
  TrackSeed_v1(const TrackSeed&);
  TrackSeed_v1(const TrackSeed_v1&);
  TrackSeed_v1& operator=(const TrackSeed_v1& seed);

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrackSeed_v1(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed&) override;
  void CopyFrom(TrackSeed* seed) override { CopyFrom(*seed); }
  PHObject* CloneMe() const override { return new TrackSeed_v1(*this); }


  ///@name accessors
  //@{
  int get_charge() const override;
  float get_pz() const override;
  float get_pt() const override;
  float get_p() const override;

  float get_eta() const override;
  float get_theta() const override;

  float get_qOverR() const override { return m_qOverR; }
  float get_X0() const override { return m_X0; }
  float get_Y0() const override { return m_Y0; }
  float get_Z0() const override { return m_Z0; }
  float get_slope() const override { return m_slope; }
  short int get_crossing() const override { return m_crossing; }

  bool empty_cluster_keys() const override { return m_cluster_keys.empty(); }

  size_t size_cluster_keys() const override { return m_cluster_keys.size(); }

  ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const override { return m_cluster_keys.find(clusterid); }
  ConstClusterKeyIter begin_cluster_keys() const override { return m_cluster_keys.begin(); }
  ConstClusterKeyIter end_cluster_keys() const override { return m_cluster_keys.end(); }
  ClusterKeyIter find_cluster_keys(unsigned int clusterid) override { return m_cluster_keys.find(clusterid); }
  ClusterKeyIter begin_cluster_keys() override { return m_cluster_keys.begin(); }
  ClusterKeyIter end_cluster_keys() override { return m_cluster_keys.end(); }
  //@}

  ///@name modifiers
  //@{
  void set_crossing(const short int crossing) override { m_crossing = crossing; }
  void set_qOverR(const float qOverR) override { m_qOverR = qOverR; }
  void set_X0(const float X0) override { m_X0 = X0; }
  void set_Y0(const float Y0) override { m_Y0 = Y0; }
  void set_Z0(const float Z0) override { m_Z0 = Z0; }
  void set_slope(const float slope) override { m_slope = slope; }

  void clear_cluster_keys() override { m_cluster_keys.clear(); }
  void insert_cluster_key(TrkrDefs::cluskey clusterid) override { m_cluster_keys.insert(clusterid); }
  size_t erase_cluster_key(TrkrDefs::cluskey clusterid) override { return m_cluster_keys.erase(clusterid); }
  //@}

 private:
  ClusterKeySet m_cluster_keys;

  float m_qOverR = NAN;
  float m_X0 = NAN;
  float m_Y0 = NAN;
  float m_slope = NAN;
  float m_Z0 = NAN;

  short int m_crossing = std::numeric_limits<short int>::max();

  ClassDefOverride(TrackSeed_v1, 1);
};

#endif
