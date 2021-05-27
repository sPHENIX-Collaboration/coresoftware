#ifndef TRACKBASEHISTORIC_SVTXTRACK_H
#define TRACKBASEHISTORIC_SVTXTRACK_H

#include "SvtxTrackState.h"

#include <trackbase/TrkrDefs.h>

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>

#include <limits.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

class SvtxTrack : public PHObject
{
 public:
  typedef std::map<float, SvtxTrackState*> StateMap;
  typedef StateMap::const_iterator ConstStateIter;
  typedef StateMap::iterator StateIter;

  typedef std::set<unsigned int> ClusterSet;
  typedef ClusterSet::const_iterator ConstClusterIter;
  typedef ClusterSet::iterator ClusterIter;

  typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
  typedef ClusterKeySet::const_iterator ConstClusterKeyIter;
  typedef ClusterKeySet::iterator ClusterKeyIter;

  enum CAL_LAYER
  {
    PRES = 0,
    CEMC = 1,
    HCALIN = 2,
    HCALOUT = 3
  };

  ~SvtxTrack() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrack base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const SvtxTrack& ) 
  {}

  //! copy content from base class
  virtual void CopyFrom( SvtxTrack* ) 
  {}

  //
  // basic track information ---------------------------------------------------
  //

  virtual unsigned int get_id() const { return UINT_MAX; }
  virtual void set_id(unsigned int id) {}

  virtual unsigned int get_vertex_id() const { return UINT_MAX; }
  virtual void set_vertex_id(unsigned int vertex_id) {}

  virtual bool get_positive_charge() const { return false; }
  virtual void set_positive_charge(bool ispos) {}

  virtual int get_charge() const { return -1; }
  virtual void set_charge(int charge) {}

  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float chisq) {}

  virtual unsigned int get_ndf() const { return UINT_MAX; }
  virtual void set_ndf(int ndf) {}

  virtual float get_quality() const { return NAN; }

  virtual float get_dca() const { return NAN; }
  virtual void set_dca(float dca) {}

  virtual float get_dca_error() const { return NAN; }
  virtual void set_dca_error(float dca) {}

  virtual float get_dca2d() const { return NAN; }
  virtual void set_dca2d(float dca2d) {}

  virtual float get_dca2d_error() const { return NAN; }
  virtual void set_dca2d_error(float error) {}

  virtual float get_dca3d_xy() const { return NAN; }
  virtual void set_dca3d_xy(float dcaxy) {}

  virtual float get_dca3d_xy_error() const { return NAN; }
  virtual void set_dca3d_xy_error(float error) {}

  virtual float get_dca3d_z() const { return NAN; }
  virtual void set_dca3d_z(float dcaz) {}

  virtual float get_dca3d_z_error() const { return NAN; }
  virtual void set_dca3d_z_error(float error) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float z) {}

  virtual float get_pos(unsigned int i) const { return NAN; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float px) {}

  virtual float get_py() const { return NAN; }
  virtual void set_py(float py) {}

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float pz) {}

  virtual float get_mom(unsigned int i) const { return NAN; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }

  virtual float get_error(int i, int j) const { return NAN; }
  virtual void set_error(int i, int j, float value) {}

  //
  // state methods -------------------------------------------------------------
  //
  virtual bool empty_states() const { return false; }
  virtual size_t size_states() const { return 0; }
  virtual size_t count_states(float pathlength) const { return 0; }
  virtual void clear_states() {}

  virtual const SvtxTrackState* get_state(float pathlength) const { return nullptr; }
  virtual SvtxTrackState* get_state(float pathlength) { return nullptr; }
  virtual SvtxTrackState* insert_state(const SvtxTrackState* state) { return nullptr; }
  virtual size_t erase_state(float pathlength) { return 0; }

  virtual ConstStateIter begin_states() const;
  virtual ConstStateIter find_state(float pathlength)  const;
  virtual ConstStateIter end_states() const;

  virtual StateIter begin_states();
  virtual StateIter find_state(float pathlength);
  virtual StateIter end_states();

  //
  // associated cluster ids methods --------------------------------------------
  //

  // needed by old tracking

  //! deprecated - please use cluster keys instead
  virtual void clear_clusters() {}
  //! deprecated - please use cluster keys instead
  virtual bool empty_clusters() const { return false; }
  //! deprecated - please use cluster keys instead
  virtual size_t size_clusters() const { return 0; }

  //! deprecated - please use cluster keys instead
  virtual void insert_cluster(unsigned int clusterid) {}
  //! deprecated - please use cluster keys instead
  virtual size_t erase_cluster(unsigned int clusterid) { return 0; }
  //! deprecated - please use cluster keys instead
  virtual ConstClusterIter begin_clusters() const;
  //! deprecated - please use cluster keys instead
  virtual ConstClusterIter find_cluster(unsigned int clusterid) const;
  //! deprecated - please use cluster keys instead
  virtual ConstClusterIter end_clusters() const;
  //! deprecated - please use cluster keys instead
  virtual ClusterIter begin_clusters();
  //! deprecated - please use cluster keys instead
  virtual ClusterIter find_cluster(unsigned int clusterid);
  //! deprecated - please use cluster keys instead
  virtual ClusterIter end_clusters();

  // needed by new tracking
  virtual void clear_cluster_keys() {}
  virtual bool empty_cluster_keys() const { return false; }
  virtual size_t size_cluster_keys() const { return 0; }

  virtual void insert_cluster_key(TrkrDefs::cluskey clusterid) {}
  virtual size_t erase_cluster_key(TrkrDefs::cluskey clusterid) { return 0; }
  virtual ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const;
  virtual ConstClusterKeyIter begin_cluster_keys() const;
  virtual ConstClusterKeyIter end_cluster_keys() const;
  virtual ClusterKeyIter begin_cluster_keys();
  virtual ClusterKeyIter find_cluster_keys(unsigned int clusterid);
  virtual ClusterKeyIter end_cluster_keys();

  //
  // calo projection methods ---------------------------------------------------
  //
  virtual float get_cal_dphi(CAL_LAYER layer) const { return 0.; }
  virtual void set_cal_dphi(CAL_LAYER layer, float dphi) {}

  virtual float get_cal_deta(CAL_LAYER layer) const { return 0.; }
  virtual void set_cal_deta(CAL_LAYER layer, float deta) {}

  virtual float get_cal_energy_3x3(CAL_LAYER layer) const { return 0.; }
  virtual void set_cal_energy_3x3(CAL_LAYER layer, float energy_3x3) {}

  virtual float get_cal_energy_5x5(CAL_LAYER layer) const { return 0.; }
  virtual void set_cal_energy_5x5(CAL_LAYER layer, float energy_5x5) {}

  virtual unsigned int get_cal_cluster_id(CAL_LAYER layer) const { return 0; }
  virtual void set_cal_cluster_id(CAL_LAYER layer, unsigned int id) {}

  virtual TrkrDefs::cluskey get_cal_cluster_key(CAL_LAYER layer) const { return 0; }
  virtual void set_cal_cluster_key(CAL_LAYER layer, TrkrDefs::cluskey key) {}

  virtual float get_cal_cluster_e(CAL_LAYER layer) const { return 0.; }
  virtual void set_cal_cluster_e(CAL_LAYER layer, float e) {}

  // Acts methods for use by Acts modules only
  virtual float get_acts_covariance(unsigned int i, unsigned int j) const { return NAN;}
  virtual void set_acts_covariance(unsigned int i, unsigned int j, float value) {}
 
  

  //
  // truth track interface ---------------------------------------------------
  //

  //SvtxTrack_FastSim
  virtual unsigned int get_truth_track_id() const { return UINT_MAX; }
  virtual void set_truth_track_id(unsigned int truthTrackId) {}
  virtual void set_num_measurements(int nmeas) {}
  virtual unsigned int get_num_measurements() const { return 0; }

  //SvtxTrack_FastSim_v1
  typedef std::map<int, std::set<PHG4HitDefs::keytype> > HitIdMap;
  typedef HitIdMap::iterator HitIdIter;
  typedef HitIdMap::const_iterator HitIdConstIter;

  virtual bool empty_g4hit_id() const { return true; }
  virtual size_t size_g4hit_id() const { return 0; }
  virtual void add_g4hit_id(int volume, PHG4HitDefs::keytype id) {}
  virtual HitIdIter begin_g4hit_id();
  virtual HitIdConstIter begin_g4hit_id() const;
  virtual HitIdIter find_g4hit_id(int volume);
  virtual HitIdConstIter find_g4hit_id(int volume) const;
  virtual HitIdIter end_g4hit_id();
  virtual HitIdConstIter end_g4hit_id() const;
  virtual size_t remove_g4hit_id(int volume, PHG4HitDefs::keytype id) { return 0; }
  virtual size_t remove_g4hit_volume(int volume) { return 0; }
  virtual void clear_g4hit_id() {}
  virtual const HitIdMap& g4hit_ids() const;

 protected:
  SvtxTrack() {}

  ClassDefOverride(SvtxTrack, 1);
};

#endif
