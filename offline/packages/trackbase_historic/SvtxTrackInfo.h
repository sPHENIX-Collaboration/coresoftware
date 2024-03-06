#ifndef TRACKBASEHISTORIC_SVTXTRACKINFO_H
#define TRACKBASEHISTORIC_SVTXTRACKINFO_H

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <limits.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>

class SvtxTrackInfo : public PHObject
{
 public:
  ~SvtxTrackInfo() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackInfo base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  virtual void CopyFrom(const SvtxTrackInfo&)
  {
  }

  //! copy content from base class
  virtual void CopyFrom(SvtxTrackInfo*)
  {
  }

  //
  // basic track information ---------------------------------------------------
  //

  virtual unsigned int get_track_id() const { return std::numeric_limits<unsigned int>::quiet_NaN(); }
  virtual void set_track_id(unsigned int ) { }

  virtual uint16_t get_subsurfkey() const { return std::numeric_limits<uint8_t>::quiet_NaN(); }
  virtual void set_subsurfkey(uint16_t ) { }

  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float) {}

  virtual uint8_t get_ndf() const { return std::numeric_limits<uint8_t>::quiet_NaN(); }
  virtual void set_ndf(uint8_t) {}

  // bitmap for clusters
  virtual uint64_t get_hitbitmap() const { return std::numeric_limits<uint64_t>::quiet_NaN(); }
  virtual void set_hitbitmap(uint64_t) {}

  // vertex information
  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_pos(unsigned int) const { return NAN; }

  virtual void set_phi(const float) {}
  virtual void set_theta(const float) {}
  virtual void set_qOp(const float) {}

  virtual float get_px() const { return NAN; }
  virtual float get_py() const { return NAN; }
  virtual float get_pz() const { return NAN; }
  virtual float get_mom(unsigned int) const { return NAN; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }
  virtual float get_theta() const { return NAN; }
  virtual float get_qOp() const { return NAN; }
  virtual int get_charge() const { return std::numeric_limits<int>::quiet_NaN(); }
  virtual float get_covariance(int /*i*/, int /*j*/) const { return NAN; }
  virtual void set_covariance(int /*i*/, int /*j*/, float /*value*/) {}

  virtual short int get_crossing() const { return SHRT_MAX; }
  virtual void set_crossing(short int /*crossing*/) {}

  virtual float get_x_outer_tpc() const { return NAN; }
  virtual void set_x_outer_tpc(float ) { }

  virtual float get_y_outer_tpc() const { return NAN; }
  virtual void set_y_outer_tpc(float ) { }

  virtual float get_z_outer_tpc() const { return NAN; }
  virtual void set_z_outer_tpc(float ) { }

  virtual void set_phi_outer_tpc(const float ) { }
  virtual void set_theta_outer_tpc(const float ) { }
  virtual void set_qOp_outer_tpc(const float ) { }

  virtual float get_px_outer_tpc() const { return NAN; }
  virtual float get_py_outer_tpc() const { return NAN; }
  virtual float get_pz_outer_tpc() const { return NAN; }

  virtual float get_mom_outer_tpc(unsigned int ) const { return NAN; }

  virtual float get_p_outer_tpc() const { return NAN; }
  virtual float get_pt_outer_tpc() const { return NAN; }
  virtual float get_eta_outer_tpc() const { return NAN; }
  virtual float get_phi_outer_tpc() const { return NAN; }
  virtual float get_theta_outer_tpc() const { return NAN; }
  virtual float get_qOp_outer_tpc() const { return NAN; }

  virtual float get_covariance_outer_tpc(int /*i*/, int /*j*/) const { return NAN;}
  virtual void set_covariance_outer_tpc(int /*i*/, int /*j*/, float /*value*/) { }

  virtual float get_x(int /*state*/) const { return NAN; }
  virtual void set_x(int /*state*/, float ) {  }

  virtual float get_y(int /*state*/) const { return NAN; }
  virtual void set_y(int /*state*/, float ) { }

  virtual float get_z(int /*state*/) const { return NAN; }
  virtual void set_z(int /*state*/, float ) { }

  virtual float get_pos(int /*state*/, unsigned int ) const { return NAN; }

  virtual float get_px(int /*state*/) const { return NAN; }
  virtual float get_py(int /*state*/) const { return NAN; }
  virtual float get_pz(int /*state*/) const { return NAN; }

  virtual void set_phi(int /*state*/, const float ) { NAN; }
  virtual void set_theta(int /*state*/, const float ) { NAN; }
  virtual void set_qOp(int /*state*/, const float ) { NAN; }

  virtual float get_mom(int /*state*/, unsigned int ) const { return NAN; }

  virtual float get_p(int /*state*/) const { return NAN; }
  virtual float get_pt(int /*state*/) const { return NAN; }
  virtual float get_eta(int /*state*/) const { return NAN; }
  virtual float get_phi(int /*state*/) const { return NAN; }
  virtual float get_theta(int /*state*/) const { return NAN; }
  virtual float get_qOp(int /*state*/) const { return NAN; }

  virtual float get_projected_eta(int /*state*/) const { return NAN; }
  virtual float get_projected_phi(int /*state*/) const { return NAN; }

  virtual float get_covariance(int /*state*/, int /*i*/, int /*j*/) const { return NAN; }
  virtual void set_covariance(int /*state*/, int /*i*/, int /*j*/, float /*value*/) { }

 private:
  ClassDefOverride(SvtxTrackInfo, 1);
};

#endif
