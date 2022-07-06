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
  TrackSeed_v1();
  
  /// Copy constructors
  TrackSeed_v1( const TrackSeed& );
  TrackSeed_v1( const TrackSeed_v1& );
  TrackSeed_v1& operator=(const TrackSeed_v1& seed);
  ~TrackSeed_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TrackSeed_v1(); }
  int isValid() const override { return 1; }
  void CopyFrom( const TrackSeed&) override;
  void CopyFrom( TrackSeed* seed) override { CopyFrom( *seed ); }
  PHObject* CloneMe() const override { return new TrackSeed_v1(*this); }

  int get_charge() const override;
  float get_px(TrkrClusterContainer *clusters,
	       ActsGeometry *tGeometry) const override;
  float get_py(TrkrClusterContainer *clusters,
	       ActsGeometry *tGeometry) const override;
  float get_pz() const override;
  float get_x() const override;
  float get_y() const override;
  float get_z() const override;
  float get_phi(TrkrClusterContainer *clusters,
		ActsGeometry *tGeometry) const override;
  float get_phi(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions) const override;
  float get_eta() const override;
  float get_theta() const override;
  float get_pt() const override;
  float get_p() const override;
  
  float get_qOverR() const override { return m_qOverR; }
  float get_X0() const override { return m_X0; }
  float get_Y0() const override { return m_Y0; }
  float get_slope() const override { return m_slope; }
  float get_Z0() const override { return m_Z0; }
  short int get_crossing() const override { return m_crossing; }

  void set_crossing(const short int crossing) override { m_crossing = crossing; }
  void set_qOverR(const float qOverR) override { m_qOverR = qOverR; }
  void set_X0(const float X0) override { m_X0 = X0; }
  void set_Y0(const float Y0) override { m_Y0 = Y0; }
  void set_slope(const float slope) override { m_slope = slope; }
  void set_Z0(const float Z0) override { m_Z0 = Z0; }

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
  void circleFitByTaubin(TrkrClusterContainer *clusters,
			 ActsGeometry *tGeometry,
			 uint8_t startLayer = 0,
			 uint8_t endLayer = 58) override;
  /// Updates r-z slope and intercept B
  void lineFit(TrkrClusterContainer *clusters,
	       ActsGeometry *tGeometry,
	       uint8_t startLayer = 0,
	       uint8_t endLayer = 58) override;
  
  void circleFitByTaubin(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
			 uint8_t startLayer = 0,
			 uint8_t endLayer = 58) override;

  void lineFit(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
	       uint8_t startLayer = 0,
	       uint8_t endLayer = 58) override;

 protected:

  /// Returns transverse PCA to (0,0)
  std::pair<float,float> findRoot() const;

 private:

  ClusterKeySet m_cluster_keys;
  
  float m_qOverR = NAN;
  float m_X0 = NAN;
  float m_Y0 = NAN;
  float m_slope = NAN;
  float m_Z0 = NAN;

  short int m_crossing = std::numeric_limits<short int>::max();

  ClassDefOverride(TrackSeed_v1,1);

};

#endif
