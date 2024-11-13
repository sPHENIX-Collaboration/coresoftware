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
  TrackSeed_v2() = default;

  /// Copy constructors
  TrackSeed_v2(const TrackSeed&);
  TrackSeed_v2(const TrackSeed_v2&);
  TrackSeed_v2& operator=(const TrackSeed_v2& seed);

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrackSeed_v2(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed&) override;
  void CopyFrom(TrackSeed* seed) override { CopyFrom(*seed); }
  PHObject* CloneMe() const override { return new TrackSeed_v2(*this); }


  ///@name accessors
  //@{
  float get_px() const override;
  float get_py() const override;
  float get_pz() const override;
  float get_p() const override;
  float get_pt() const override;

  float get_eta() const override;
  float get_theta() const override;


  //methods that return member variables
  int get_charge() const override;
  float get_qOverR() const override { return m_qOverR; }
  float get_X0() const override { return m_X0; }
  float get_Y0() const override { return m_Y0; }
  float get_Z0() const override { return m_Z0; }
  float get_slope() const override { return m_slope; }
  float get_phi() const override  { return m_phi; }  // returns the stored phi
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

  ///@modifiers
  //@{

  void set_crossing(const short int crossing) override { m_crossing = crossing; }
  void set_qOverR(const float qOverR) override { m_qOverR = qOverR; }
  void set_X0(const float X0) override { m_X0 = X0; }
  void set_Y0(const float Y0) override { m_Y0 = Y0; }
  void set_Z0(const float Z0) override { m_Z0 = Z0; }
  void set_slope(const float slope) override { m_slope = slope; }
  void set_phi(const float phi) override { m_phi = phi; }

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
  float m_phi = NAN;

  short int m_crossing = std::numeric_limits<short int>::max();

  ClassDefOverride(TrackSeed_v2, 1);
};

#endif
