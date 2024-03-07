#ifndef TRACKBASEHISTORIC_SVTXTRACKINFOV3_H
#define TRACKBASEHISTORIC_SVTXTRACKINFOV3_H

#include "SvtxTrackInfo.h"
#include "TrackStateInfo_v1.h"

class PHObject;

class SvtxTrackInfo_v3: public SvtxTrackInfo
{
 public:
  enum STATE
  {
    VERTEX = 0,
    OUTER_TPC = 1,
    EMCAL_FRONTFACE = 2,
    EMCAL_BACKFACE = 3,
    HCALIN_FRONTFACE = 4,
    HCALIN_BACKFACE = 5,
    HCALOUT_FRONTFACE = 6,
    HCALOUT_BACKFACE = 7
  };

  SvtxTrackInfo_v3() : m_states(8, TrackStateInfo_v1()) {}

  //* base class copy constructor
  SvtxTrackInfo_v3( const SvtxTrackInfo& ) : m_states(8, TrackStateInfo_v1()) {}

  //* copy constructor
  SvtxTrackInfo_v3(const SvtxTrackInfo_v3& source) : m_states(8, TrackStateInfo_v1())
  {
    m_track_id = source.get_track_id();
    m_outer_tpc_subsurfkey = source.get_subsurfkey();
    m_chisq = source.get_chisq();
    m_ndf = source.get_ndf();
    m_crossing = source.get_crossing();
    m_hitbitmap = source.get_hitbitmap();

    for(int istate = STATE::VERTEX; istate <= STATE::HCALOUT_BACKFACE; istate++)
    {
      set_x(istate, source.get_x(istate));
      set_y(istate, source.get_y(istate));
      set_z(istate, source.get_z(istate));
      set_phi(istate, source.get_phi());
      set_theta(istate, source.get_theta());
      set_qOp(istate, source.get_qOp());

      for (int i = 0; i < 5; i++)
      {
        for (int j = i; j < 5; j++)
        {
          set_covariance(istate, i, j, source.get_covariance(istate, i, j));
        }
      }
    }
  }

  //* assignment operator
  SvtxTrackInfo_v3& operator=(const SvtxTrackInfo_v3& track);

  //* destructor
  ~SvtxTrackInfo_v3() override {}

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackInfo_v3 class" << std::endl;
  }
  void Reset() override { *this = SvtxTrackInfo_v3(); }

  PHObject* CloneMe() const override { return new SvtxTrackInfo_v3(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content
  void CopyFrom( const SvtxTrackInfo& ) override;
  void CopyFrom( SvtxTrackInfo* source ) override
  {
    CopyFrom( *source );
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

  float get_x() const override { return m_states.at(STATE::VERTEX).get_x(); }
  void set_x(float x) override { m_states.at(STATE::VERTEX).set_x(x); }

  float get_y() const override { return m_states.at(STATE::VERTEX).get_y(); }
  void set_y(float y) override { m_states.at(STATE::VERTEX).set_y(y); }

  float get_z() const override { return m_states.at(STATE::VERTEX).get_z(); }
  void set_z(float z) override { m_states.at(STATE::VERTEX).set_z(z); }

  float get_pos(unsigned int i) const override { return m_states.at(STATE::VERTEX).get_pos(i); }

  float get_px() const override { return m_states.at(STATE::VERTEX).get_px(); }
  float get_py() const override { return m_states.at(STATE::VERTEX).get_py(); }
  float get_pz() const override { return m_states.at(STATE::VERTEX).get_pz(); }

  void set_phi(const float phi) override { m_states.at(STATE::VERTEX).set_phi(phi); }
  void set_theta(const float theta) override { m_states.at(STATE::VERTEX).set_theta(theta); }
  void set_qOp(const float qop) override { m_states.at(STATE::VERTEX).set_qOp(qop); }

  float get_mom(unsigned int i) const override { return m_states.at(STATE::VERTEX).get_mom(i); }

  float get_p() const override { return m_states.at(STATE::VERTEX).get_p(); }
  float get_pt() const override { return m_states.at(STATE::VERTEX).get_pt(); }
  float get_eta() const override { return m_states.at(STATE::VERTEX).get_eta(); }
  float get_phi() const override { return m_states.at(STATE::VERTEX).get_phi(); }
  float get_theta() const override { return m_states.at(STATE::VERTEX).get_theta(); }
  float get_qOp() const override { return m_states.at(STATE::VERTEX).get_qOp(); }
  int get_charge() const override { return m_states.at(STATE::VERTEX).get_charge(); }

  float get_covariance(int i, int j) const override { return m_states.at(STATE::VERTEX).get_covariance(i, j);}
  void set_covariance(int i, int j, float value) override {m_states.at(STATE::VERTEX).set_covariance(i, j, value);}

  float get_x(int state) const override { return m_states.at(state).get_x(); }
  void set_x(int state, float x) override { m_states.at(state).set_x(x); }

  float get_y(int state) const override { return m_states.at(state).get_y(); }
  void set_y(int state, float y) override { m_states.at(state).set_y(y); }

  float get_z(int state) const override { return m_states.at(state).get_z(); }
  void set_z(int state, float z) override { m_states.at(state).set_z(z); }

  float get_pos(int state, unsigned int i) const override { return m_states.at(state).get_pos(i); }

  float get_px(int state) const override { return m_states.at(state).get_px(); }
  float get_py(int state) const override { return m_states.at(state).get_py(); }
  float get_pz(int state) const override { return m_states.at(state).get_pz(); }

  void set_phi(int state, const float phi) override { m_states.at(state).set_phi(phi); }
  void set_theta(int state, const float theta) override { m_states.at(state).set_theta(theta); }
  void set_qOp(int state, const float qop) override { m_states.at(state).set_qOp(qop); }

  float get_mom(int state, unsigned int i) const override { return m_states.at(state).get_mom(i); }

  float get_p(int state) const override { return m_states.at(state).get_p(); }
  float get_pt(int state) const override { return m_states.at(state).get_pt(); }
  float get_eta(int state) const override { return m_states.at(state).get_eta(); }
  float get_phi(int state) const override { return m_states.at(state).get_phi(); }
  float get_theta(int state) const override { return m_states.at(state).get_theta(); }
  float get_qOp(int state) const override { return m_states.at(state).get_qOp(); }

  float get_projected_eta(int state) const override { return asinh(m_states.at(state).get_z()/sqrt(pow(m_states.at(state).get_x(), 2) + pow(m_states.at(state).get_y(), 2))); }
  float get_projected_phi(int state) const override { return atan2(m_states.at(state).get_y(), m_states.at(state).get_x()); }

  float get_covariance(int state, int i, int j) const override { return m_states.at(state).get_covariance(i, j);}
  void set_covariance(int state, int i, int j, float value) override {m_states.at(state).set_covariance(i, j, value);}

 private:

  // track information
  unsigned int m_track_id = std::numeric_limits<unsigned int>::quiet_NaN();
  uint16_t m_outer_tpc_subsurfkey = std::numeric_limits<uint16_t>::quiet_NaN();
  float m_chisq = std::numeric_limits<float>::quiet_NaN();
  uint8_t m_ndf = std::numeric_limits<uint8_t>::quiet_NaN();
  uint64_t m_hitbitmap = std::numeric_limits<uint64_t>::quiet_NaN();
  short int m_crossing = std::numeric_limits<short int>::quiet_NaN();

  // track state information
  std::vector<TrackStateInfo_v1> m_states;

  ClassDefOverride(SvtxTrackInfo_v3, 1)
};

#endif
