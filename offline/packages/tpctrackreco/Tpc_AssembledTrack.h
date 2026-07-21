#pragma once

#include <trackbase/TrkrDefs.h>
#include <phool/PHObject.h>

#include <iostream>
#include <utility>
#include <vector>

// A full TPC track built by connecting already reconstructed Tpc_ModuleTrack
// pieces across modules/sectors.  Concrete payload lives in versioned
// subclasses (Tpc_AssembledTrackv1, ...), following the same schema-evolution pattern
// as Tpc_ModuleTrack/Tpc_ModuleTrackv1.
class Tpc_AssembledTrack : public PHObject
{
 public:
  Tpc_AssembledTrack() = default;
  ~Tpc_AssembledTrack() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_AssembledTrack base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  // --- Identity ---
  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual int get_side() const { return 0; }

  // --- Topology ---
  virtual unsigned int get_nsegments() const { return 0; }
  virtual unsigned int get_nblobs() const { return 0; }
  virtual unsigned int get_nrawhits() const { return 0; }
  virtual unsigned int get_first_layer() const { return 0; }
  virtual unsigned int get_last_layer() const { return 0; }
  virtual unsigned int get_first_sector() const { return 0; }
  virtual unsigned int get_last_sector() const { return 0; }
  virtual unsigned int get_first_region() const { return 0; }
  virtual unsigned int get_last_region() const { return 0; }

  // --- Full-track fit in global sector-unwrapped coordinates ---
  // global_phi is in radians.  It is built from sector + local pad using the
  // same local pad/phi calibration constants as Tpc_ModuleTrackReco.
  virtual double get_phi_slope() const { return 0.0; }
  virtual double get_phi_intercept() const { return 0.0; }
  virtual double get_tbin_slope() const { return 0.0; }
  virtual double get_tbin_intercept() const { return 0.0; }
  virtual double get_chi2_phi() const { return 0.0; }
  virtual double get_chi2_tbin() const { return 0.0; }
  virtual int get_ndof_phi() const { return 0; }
  virtual int get_ndof_tbin() const { return 0; }

  // --- Setters ---
  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_side(int) {}
  virtual void set_nsegments(unsigned int) {}
  virtual void set_nblobs(unsigned int) {}
  virtual void set_nrawhits(unsigned int) {}
  virtual void set_first_layer(unsigned int) {}
  virtual void set_last_layer(unsigned int) {}
  virtual void set_first_sector(unsigned int) {}
  virtual void set_last_sector(unsigned int) {}
  virtual void set_first_region(unsigned int) {}
  virtual void set_last_region(unsigned int) {}
  virtual void set_phi_slope(double) {}
  virtual void set_phi_intercept(double) {}
  virtual void set_tbin_slope(double) {}
  virtual void set_tbin_intercept(double) {}
  virtual void set_chi2_phi(double) {}
  virtual void set_chi2_tbin(double) {}
  virtual void set_ndof_phi(int) {}
  virtual void set_ndof_tbin(int) {}

  // --- Event vertex estimated from full-track pairs ---
  // r/phi are from the r-phi projection; tbin is from the r-tbin projection.
  virtual int get_vertex_valid() const { return 0; }
  virtual double get_vertex_x() const { return 0.0; }
  virtual double get_vertex_y() const { return 0.0; }
  virtual double get_vertex_r() const { return 0.0; }
  virtual double get_vertex_phi() const { return 0.0; }
  virtual double get_vertex_tbin() const { return 0.0; }
  virtual unsigned int get_vertex_npairs() const { return 0; }
  virtual double get_vertex_quality() const { return 0.0; }

  virtual void set_vertex_valid(int) {}
  virtual void set_vertex_x(double) {}
  virtual void set_vertex_y(double) {}
  virtual void set_vertex_r(double) {}
  virtual void set_vertex_phi(double) {}
  virtual void set_vertex_tbin(double) {}
  virtual void set_vertex_npairs(unsigned int) {}
  virtual void set_vertex_quality(double) {}

  // --- Preliminary GenFit seed from the connector fit ---
  virtual int get_seed_valid() const { return 0; }
  virtual double get_seed_x() const { return 0.0; }
  virtual double get_seed_y() const { return 0.0; }
  virtual double get_seed_z() const { return 0.0; }
  virtual double get_seed_px() const { return 0.0; }
  virtual double get_seed_py() const { return 0.0; }
  virtual double get_seed_pz() const { return 0.0; }
  virtual double get_seed_cov(unsigned int, unsigned int) const { return 0.0; }

  virtual void set_seed_valid(int) {}
  virtual void set_seed_x(double) {}
  virtual void set_seed_y(double) {}
  virtual void set_seed_z(double) {}
  virtual void set_seed_px(double) {}
  virtual void set_seed_py(double) {}
  virtual void set_seed_pz(double) {}
  virtual void set_seed_cov(unsigned int, unsigned int, double) {}
  
  // Source Tpc_ModuleTrack ids used by this assembled track.
  virtual void add_source_track(unsigned int /*track_id*/,
                                unsigned int /*region*/,
                                unsigned int /*sector*/) {}
  virtual unsigned int size_source_tracks() const { return 0; }
  virtual unsigned int get_source_track_id(unsigned int) const { return 0; }
  virtual unsigned int get_source_region(unsigned int) const { return 0; }
  virtual unsigned int get_source_sector(unsigned int) const { return 0; }

  // Hit index references copied from source in-module tracks.
  using HitIndex = std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>;
  virtual void add_hit_index(TrkrDefs::hitsetkey, TrkrDefs::hitkey) {}
  virtual unsigned int size_hit_indices() const { return 0; }
  virtual HitIndex get_hit_index(unsigned int) const { return {0, 0}; }
  virtual const std::vector<HitIndex>& get_hit_indices() const
  {
    static const std::vector<HitIndex> empty_indices;
    return empty_indices;
  }

 private:
  ClassDefOverride(Tpc_AssembledTrack, 0)
};
