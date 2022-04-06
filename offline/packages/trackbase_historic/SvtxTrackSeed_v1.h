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

  int get_charge() const override { return (get_positive_charge()) ? 1 : -1; }
  bool get_positive_charge() const override { return m_is_positive_charge; }
  float get_px() const override { return m_px; }
  float get_py() const override { return m_py; }
  float get_pz() const override { return m_pz; }
  float get_x() const override { return m_pcax; }
  float get_y() const override { return m_pcay; }
  float get_z() const override { return m_pcaz; }
  float get_phi() const override { return atan2(m_py,m_px); }
  float get_eta() const override { return asinh(m_pz/get_pt()); }
  float get_pt() const override { return sqrt(m_px*m_px + m_py*m_py); }
  float get_p() const override { return sqrt(m_px*m_px + m_py*m_py + m_pz*m_pz); }
  
  void set_positive_charge(bool pos) override { m_is_positive_charge = pos; }
  void set_charge(const int charge) override { (charge > 0) ? set_positive_charge(true) : set_positive_charge(false); }
  void set_px(const float px) override { m_px = px; }
  void set_py(const float py) override { m_py = py; }
  void set_pz(const float pz) override { m_pz = pz; }
  void set_x(const float pcax) override { m_pcax = pcax; }
  void set_y(const float pcay) override { m_pcay = pcay; }
  void set_z(const float pcaz) override { m_pcaz = pcaz; }

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


 private:
  
  bool m_is_positive_charge = false;
  ClusterKeySet m_cluster_keys;
  
  float m_px = NAN;
  float m_py = NAN;
  float m_pz = NAN;
  float m_pcax = NAN;
  float m_pcay = NAN;
  float m_pcaz = NAN;

  ClassDefOverride(SvtxTrackSeed_v1,1);

};

#endif
