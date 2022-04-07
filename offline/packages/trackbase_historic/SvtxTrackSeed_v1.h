#ifndef TRACKBASEHISTORIC_SVTXTRACKSEEDV1_H
#define TRACKBASEHISTORIC_SVTXTRACKSEEDV1_H

#include "SvtxTrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <limits.h>
#include <cmath>
#include <iostream>

class SvtxTrackSeed_v1 : public SvtxTrackSeed
{
 public: 
  SvtxTrackSeed_v1();
  
  /// Copy constructors
  SvtxTrackSeed_v1( const SvtxTrackSeed& );
  SvtxTrackSeed_v1( const SvtxTrackSeed_v1& );
  SvtxTrackSeed_v1& operator=(const SvtxTrackSeed_v1& seed);
  ~SvtxTrackSeed_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrackSeed_v1(); }
  void CopyFrom( const SvtxTrackSeed&) override;
  void CopyFrom( SvtxTrackSeed* seed) override { CopyFrom( *seed ); }

  int get_charge() const override;
  float get_px() const override;
  float get_py() const override;
  float get_pz() const override;
  float get_x() const override;
  float get_y() const override;
  float get_z() const override;
  float get_phi() const override;
  float get_eta() const override;
  float get_theta() const override;
  float get_pt() const override;
  float get_p() const override;
  
  float get_qOverR() const override { return m_qOverR; }
  float get_X0() const override { return m_X0; }
  float get_Y0() const override { return m_Y0; }
  float get_slope() const override { return m_slope; }
  float get_B() const override { return m_B; }

  void set_qOverR(const float qOverR) override { m_qOverR = qOverR; }
  void set_X0(const float X0) override { m_X0 = X0; }
  void set_Y0(const float Y0) override { m_Y0 = Y0; }
  void set_slope(const float slope) override { m_slope = slope; }
  void set_B(const float B) override { m_B = B; }

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
			 ActsSurfaceMaps *surfMaps, 
			 ActsTrackingGeometry *tGeometry) override;
  /// Updates r-z slope and intercept B
  void lineFit(TrkrClusterContainer *clusters,
	       ActsSurfaceMaps *surfMaps, 
	       ActsTrackingGeometry *tGeometry) override;
  
 private:
  /// Returns transverse PCA to (0,0)
  void findRoot(float& x ,float& y) const;
  float findRoot(bool findX) const;
  ClusterKeySet m_cluster_keys;
  
  float m_qOverR = NAN;
  float m_X0 = NAN;
  float m_Y0 = NAN;
  float m_slope = NAN;
  float m_B = NAN;

  ClassDefOverride(SvtxTrackSeed_v1,1);

};

#endif
