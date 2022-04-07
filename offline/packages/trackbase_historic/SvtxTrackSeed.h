#ifndef TRACKBASEHISTORIC_SVTXTRACKSEED_H
#define TRACKBASEHISTORIC_SVTXTRACKSEED_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/TrkrClusterContainer.h>

#include <cmath>
#include <iostream>
#include <set>

class SvtxTrackSeed : public PHObject
{
 public: 
  typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
  typedef ClusterKeySet::const_iterator ConstClusterKeyIter;
  typedef ClusterKeySet::iterator ClusterKeyIter;

  ~SvtxTrackSeed() override = default;
  
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackSeed base class\n";
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  using PHObject::CopyFrom;
  
  virtual void CopyFrom( const SvtxTrackSeed& ) {}

  virtual void CopyFrom( SvtxTrackSeed* ) {}

  virtual int get_charge() const { return -99999; }
  virtual float get_px() const { return NAN; }
  virtual float get_py() const { return NAN; }
  virtual float get_pz() const { return NAN; }
  virtual float get_x() const { return NAN; }
  virtual float get_y() const { return NAN; }
  virtual float get_z() const { return NAN; }
  virtual float get_qOverR() const { return NAN; }
  virtual float get_X0() const { return NAN; }
  virtual float get_Y0() const { return NAN; }
  virtual float get_slope() const { return NAN; }
  virtual float get_B() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_theta() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_p() const { return NAN; }
  
  virtual void set_qOverR(const float) {}
  virtual void set_X0(const float) {}
  virtual void set_Y0(const float) {}
  virtual void set_slope(const float) {}
  virtual void set_B(const float) {}

  virtual void circleFitByTaubin(TrkrClusterContainer *clusters,
				 ActsSurfaceMaps* surfMaps,
				 ActsTrackingGeometry* tGeometry) {}
  virtual void lineFit(TrkrClusterContainer *clusters,
		       ActsSurfaceMaps* surfMaps,
		       ActsTrackingGeometry* tGeometry) {}
  

  virtual void clear_cluster_keys() {}
  virtual bool empty_cluster_keys() const { return false; }
  virtual size_t size_cluster_keys() const { return 0; }
  virtual void insert_cluster_key(TrkrDefs::cluskey /*clusterid*/) {}
  virtual size_t erase_cluster_key(TrkrDefs::cluskey /*clusterid*/) { return 0; }
  virtual ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const;
  virtual ConstClusterKeyIter begin_cluster_keys() const;
  virtual ConstClusterKeyIter end_cluster_keys() const;
  virtual ClusterKeyIter begin_cluster_keys();
  virtual ClusterKeyIter find_cluster_keys(unsigned int clusterid);
  virtual ClusterKeyIter end_cluster_keys();

 protected:
  SvtxTrackSeed() {}

  ClassDefOverride(SvtxTrackSeed, 1);

};

#endif 
