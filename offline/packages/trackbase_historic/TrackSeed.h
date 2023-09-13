#ifndef TRACKBASEHISTORIC_TRACKSEED_H
#define TRACKBASEHISTORIC_TRACKSEED_H

#include <phool/PHObject.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4main/PHG4HitDefs.h>

#include <cmath>
#include <iostream>
#include <set>
#include <iterator>

class TrackSeed : public PHObject
{
 public: 
  typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
  typedef ClusterKeySet::const_iterator ConstClusterKeyIter;
  typedef ClusterKeySet::iterator ClusterKeyIter;

  ~TrackSeed() override = default;
  
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackSeed base class\n";
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  using PHObject::CopyFrom;
  
  virtual void CopyFrom( const TrackSeed& ) {}

  virtual void CopyFrom( TrackSeed* ) {}

  virtual int get_charge() const { return std::numeric_limits<int>::max(); }
  /// We need access to the first two clusters to get the phi angle right
  virtual float get_px(TrkrClusterContainer*,
		       ActsGeometry*) const { return NAN; }
  virtual float get_py(TrkrClusterContainer*,
		       ActsGeometry*) const { return NAN; }
  virtual float get_phi(TrkrClusterContainer*,
			ActsGeometry*) const { return NAN; }
  virtual float get_phi(std::map<TrkrDefs::cluskey, Acts::Vector3>&) const { return NAN; }
  virtual float get_pz() const { return NAN; }
  virtual float get_x() const { return NAN; }
  virtual float get_y() const { return NAN; }
  virtual float get_z() const { return NAN; }
  virtual float get_qOverR() const { return NAN; }
  virtual float get_X0() const { return NAN; }
  virtual float get_Y0() const { return NAN; }
  virtual float get_slope() const { return NAN; }
  virtual float get_Z0() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_theta() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_p() const { return NAN; }
  virtual short int get_crossing() const { return 0; }

  virtual void set_crossing(const short int) {}
  virtual void set_qOverR(const float) {}
  virtual void set_X0(const float) {}
  virtual void set_Y0(const float) {}
  virtual void set_slope(const float) {}
  virtual void set_Z0(const float) {}

  virtual void circleFitByTaubin(TrkrClusterContainer*,
				 ActsGeometry*,
				 uint8_t, uint8_t) {}
  virtual void lineFit(TrkrClusterContainer*,
		       ActsGeometry*,
		       uint8_t, uint8_t) {}
  
  /// In case the global cluster positions have already been obtained, 
  /// these can be called to avoid performing transformations twice
  virtual void circleFitByTaubin(std::map<TrkrDefs::cluskey, Acts::Vector3>&, 
				 uint8_t, uint8_t) {}
  virtual void lineFit(std::map<TrkrDefs::cluskey, Acts::Vector3>&, 
		       uint8_t, uint8_t) {}

  virtual void clear_cluster_keys() {}
  virtual bool empty_cluster_keys() const { return false; }
  virtual size_t size_cluster_keys() const { return 0; }
  virtual void insert_cluster_key(TrkrDefs::cluskey) {}
  virtual size_t erase_cluster_key(TrkrDefs::cluskey) { return 0; }
  virtual ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey) const;
  virtual ConstClusterKeyIter begin_cluster_keys() const;
  virtual ConstClusterKeyIter end_cluster_keys() const;
  virtual ClusterKeyIter begin_cluster_keys();
  virtual ClusterKeyIter find_cluster_keys(unsigned int);
  virtual ClusterKeyIter end_cluster_keys();

  virtual void set_silicon_seed_index(const unsigned int) {}
  virtual void set_tpc_seed_index(const unsigned int) {}
  virtual unsigned int get_silicon_seed_index() const { return 0; }
  virtual unsigned int get_tpc_seed_index() const { return 0; }

  /* ---------------------------------------------------------
   * Truth tracking interfaces
   * ---------------------------------------------------------
   */
  typedef std::map<int, std::set<PHG4HitDefs::keytype> > HitIdMap;
  typedef HitIdMap::iterator HitIdIter;
  typedef HitIdMap::const_iterator HitIdConstIter;

  virtual unsigned int get_truth_track_id() const { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_truth_track_id(unsigned int) {}
  virtual void set_num_measurements(int) {}
  virtual unsigned int get_num_measurements() const { return 0; }
  virtual bool empty_g4hit_id() const { return true; }
  virtual size_t size_g4hit_id() const { return 0; }
  virtual void add_g4hit_id(int, PHG4HitDefs::keytype) {}
  virtual HitIdIter begin_g4hit_id();
  virtual HitIdConstIter begin_g4hit_id() const;
  virtual HitIdIter find_g4hit_id(int);
  virtual HitIdConstIter find_g4hit_id(int) const;
  virtual HitIdIter end_g4hit_id();
  virtual HitIdConstIter end_g4hit_id() const;
  virtual size_t remove_g4hit_id(int, PHG4HitDefs::keytype) { return 0; }
  virtual size_t remove_g4hit_volume(int) { return 0; }
  virtual void clear_g4hit_id() {}
  virtual const HitIdMap& g4hit_ids() const;


 protected:
  TrackSeed() {}

  ClassDefOverride(TrackSeed, 1);

};

#endif 
