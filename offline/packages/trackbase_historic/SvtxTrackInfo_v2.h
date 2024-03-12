#ifndef TRACKBASEHISTORIC_SVTXTRACKINFOV2_H
#define TRACKBASEHISTORIC_SVTXTRACKINFOV2_H

#include "SvtxTrackInfo.h"
#include "TrackStateInfo_v1.h"

class PHObject;

class SvtxTrackInfo_v2 final : public SvtxTrackInfo
{
 public:
  SvtxTrackInfo_v2() {}

  //* base class copy constructor
  SvtxTrackInfo_v2(const SvtxTrackInfo&) {}

  //* copy constructor
  SvtxTrackInfo_v2(const SvtxTrackInfo_v2& source)
  {
    m_track_id = source.get_track_id();
    m_outer_tpc_subsurfkey = source.get_subsurfkey();
    m_chisq = source.get_chisq();
    m_ndf = source.get_ndf();
    m_crossing = source.get_crossing();
    m_hitbitmap = source.get_hitbitmap();

    set_x(source.get_x());
    set_y(source.get_y());
    set_z(source.get_z());
    set_phi(source.get_phi());
    set_theta(source.get_theta());
    set_qOp(source.get_qOp());

    set_x_outer_tpc(source.get_x_outer_tpc());
    set_y_outer_tpc(source.get_y_outer_tpc());
    set_z_outer_tpc(source.get_z_outer_tpc());
    set_phi_outer_tpc(source.get_phi_outer_tpc());
    set_theta_outer_tpc(source.get_theta_outer_tpc());
    set_qOp_outer_tpc(source.get_qOp_outer_tpc());

    for (int i = 0; i < 5; i++)
    {
      for (int j = i; j < 5; j++)
      {
        set_covariance(i, j, source.get_covariance(i, j));
        set_covariance_outer_tpc(i, j, source.get_covariance_outer_tpc(i, j));
      }
    }
  }

  //* assignment operator
  SvtxTrackInfo_v2& operator=(const SvtxTrackInfo_v2& track);

  //* destructor
  ~SvtxTrackInfo_v2() override {}

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackInfo_v2 class" << std::endl;
  }
  void Reset() override { *this = SvtxTrackInfo_v2(); }

  PHObject* CloneMe() const override { return new SvtxTrackInfo_v2(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content
  void CopyFrom(const SvtxTrackInfo&) override;
  void CopyFrom(SvtxTrackInfo* source) override
  {
    CopyFrom(*source);
  }

  //
  // basic track information ---------------------------------------------------
  //

  unsigned int get_track_id() const override { return m_track_id; }
  void set_track_id(unsigned int id) override { m_track_id = id; }

  uint16_t get_subsurfkey() const override { return m_outer_tpc_subsurfkey; }
  void set_subsurfkey(uint16_t key) override { m_outer_tpc_subsurfkey = key; }

  float get_chisq() const override { return m_chisq; }
  void set_chisq(float chisq) override { m_chisq = chisq; }

  uint8_t get_ndf() const override { return m_ndf; }
  void set_ndf(uint8_t ndf) override { m_ndf = ndf; }

  uint64_t get_hitbitmap() const override { return m_hitbitmap; }
  void set_hitbitmap(uint64_t hitbitmap) override { m_hitbitmap = hitbitmap; }

  short int get_crossing() const override { return m_crossing; }
  void set_crossing(short int crossing) override { m_crossing = crossing; }

  float get_x() const override { return m_state_vertex.get_x(); }
  void set_x(float x) override { m_state_vertex.set_x(x); }

  float get_y() const override { return m_state_vertex.get_y(); }
  void set_y(float y) override { m_state_vertex.set_y(y); }

  float get_z() const override { return m_state_vertex.get_z(); }
  void set_z(float z) override { m_state_vertex.set_z(z); }

  void set_phi(const float phi) override { m_state_vertex.set_phi(phi); }
  void set_theta(const float theta) override { m_state_vertex.set_theta(theta); }
  void set_qOp(const float qop) override { m_state_vertex.set_qOp(qop); }

  float get_pos(unsigned int i) const override { return m_state_vertex.get_pos(i); }

  float get_px() const override { return m_state_vertex.get_px(); }
  float get_py() const override { return m_state_vertex.get_py(); }
  float get_pz() const override { return m_state_vertex.get_pz(); }
  float get_mom(unsigned int i) const override { return m_state_vertex.get_mom(i); }

  float get_p() const override { return m_state_vertex.get_p(); }
  float get_pt() const override { return m_state_vertex.get_pt(); }
  float get_eta() const override { return m_state_vertex.get_eta(); }
  float get_phi() const override { return m_state_vertex.get_phi(); }
  float get_theta() const override { return m_state_vertex.get_theta(); }
  float get_qOp() const override { return m_state_vertex.get_qOp(); }
  int get_charge() const override { return m_state_vertex.get_charge(); }

  float get_covariance(int i, int j) const override { return m_state_vertex.get_covariance(i, j); }
  void set_covariance(int i, int j, float value) override { m_state_vertex.set_covariance(i, j, value); }

  float get_x_outer_tpc() const override { return m_state_outer_tpc.get_x(); }
  void set_x_outer_tpc(float x) override { m_state_outer_tpc.set_x(x); }

  float get_y_outer_tpc() const override { return m_state_outer_tpc.get_y(); }
  void set_y_outer_tpc(float y) override { m_state_outer_tpc.set_y(y); }

  float get_z_outer_tpc() const override { return m_state_outer_tpc.get_z(); }
  void set_z_outer_tpc(float z) override { m_state_outer_tpc.set_z(z); }

  void set_phi_outer_tpc(const float phi) override { m_state_outer_tpc.set_phi(phi); }
  void set_theta_outer_tpc(const float theta) override { m_state_outer_tpc.set_theta(theta); }
  void set_qOp_outer_tpc(const float qop) override { m_state_outer_tpc.set_qOp(qop); }

  float get_px_outer_tpc() const override { return m_state_outer_tpc.get_px(); }
  float get_py_outer_tpc() const override { return m_state_outer_tpc.get_py(); }
  float get_pz_outer_tpc() const override { return m_state_outer_tpc.get_pz(); }

  float get_mom_outer_tpc(unsigned int i) const override { return m_state_outer_tpc.get_mom(i); }

  float get_p_outer_tpc() const override { return m_state_outer_tpc.get_p(); }
  float get_pt_outer_tpc() const override { return m_state_outer_tpc.get_pt(); }
  float get_eta_outer_tpc() const override { return m_state_outer_tpc.get_eta(); }
  float get_phi_outer_tpc() const override { return m_state_outer_tpc.get_phi(); }
  float get_theta_outer_tpc() const override { return m_state_outer_tpc.get_theta(); }
  float get_qOp_outer_tpc() const override { return m_state_outer_tpc.get_qOp(); }

  float get_covariance_outer_tpc(int i, int j) const override { return m_state_outer_tpc.get_covariance(i, j); }
  void set_covariance_outer_tpc(int i, int j, float value) override { m_state_outer_tpc.set_covariance(i, j, value); }

 private:
  // track information
  unsigned int m_track_id = std::numeric_limits<unsigned int>::quiet_NaN();
  uint16_t m_outer_tpc_subsurfkey = std::numeric_limits<uint16_t>::quiet_NaN();
  float m_chisq = std::numeric_limits<float>::quiet_NaN();
  uint8_t m_ndf = std::numeric_limits<uint8_t>::quiet_NaN();
  uint64_t m_hitbitmap = std::numeric_limits<uint64_t>::quiet_NaN();
  short int m_crossing = std::numeric_limits<short int>::quiet_NaN();

  // track state information
  TrackStateInfo_v1 m_state_vertex;
  TrackStateInfo_v1 m_state_outer_tpc;

  ClassDefOverride(SvtxTrackInfo_v2, 1)
};

#endif
