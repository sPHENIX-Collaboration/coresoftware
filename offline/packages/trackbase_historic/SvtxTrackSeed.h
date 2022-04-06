#ifndef TRACKBASEHISTORIC_SVTXTRACKSEED_H
#define TRACKBASEHISTORIC_SVTXTRACKSEED_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

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
  virtual bool get_positive_charge() const { return false; }
  virtual float get_px() const { return NAN; }
  virtual float get_py() const { return NAN; }
  virtual float get_pz() const { return NAN; }
  virtual float get_x() const { return NAN; }
  virtual float get_y() const { return NAN; }
  virtual float get_z() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_p() const { return NAN; }
  
  virtual void set_positive_charge(bool) {}
  virtual void set_charge(const int) {}
  virtual void set_px(const float) {}
  virtual void set_py(const float) {}
  virtual void set_pz(const float) {}
  virtual void set_x(const float) {}
  virtual void set_y(const float) {}
  virtual void set_z(const float) {}


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
