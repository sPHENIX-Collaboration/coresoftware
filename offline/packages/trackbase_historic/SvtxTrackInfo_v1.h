#ifndef TRACKBASEHISTORIC_SVTXTRACKINFOV1_H
#define TRACKBASEHISTORIC_SVTXTRACKINFOV1_H

#include "SvtxTrackInfo.h"
#include "TrackStateInfo_v1.h"

#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>  // for pair

class PHObject;

class SvtxTrackInfo_v1 : public SvtxTrackInfo
{
 public:
  SvtxTrackInfo_v1() {}

  //* base class copy constructor
  SvtxTrackInfo_v1(const SvtxTrackInfo&) {}

  //* copy constructor
  SvtxTrackInfo_v1(const SvtxTrackInfo_v1& source)
  {
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

    for (int i = 0; i < 6; i++)
    {
      for (int j = i; j < 6; j++)
      {
        set_covariance(i, j, source.get_covariance(i, j));
      }
    }
  }

  //* assignment operator
  SvtxTrackInfo_v1& operator=(const SvtxTrackInfo_v1& track);

  //* destructor
  ~SvtxTrackInfo_v1() override {}

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackInfo_v1 class" << std::endl;
  }
  void Reset() override { *this = SvtxTrackInfo_v1(); }
  // int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxTrackInfo_v1(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content from base class
  void CopyFrom(const SvtxTrackInfo&) override;
  void CopyFrom(SvtxTrackInfo* source) override
  {
    CopyFrom(*source);
  }

  //
  // basic track information ---------------------------------------------------
  //

  float get_chisq() const override { return m_chisq; }
  void set_chisq(float chisq) override { m_chisq = chisq; }

  uint8_t get_ndf() const override { return m_ndf; }
  void set_ndf(uint8_t ndf) override { m_ndf = ndf; }

  uint64_t get_hitbitmap() const override { return m_hitbitmap; }
  void set_hitbitmap(uint64_t hitbitmap) override { m_hitbitmap = hitbitmap; }

  short int get_crossing() const override { return m_crossing; }
  void set_crossing(short int crossing) override { m_crossing = crossing; }

  float get_x() const override { return m_state.get_x(); }
  void set_x(float x) override { m_state.set_x(x); }

  float get_y() const override { return m_state.get_y(); }
  void set_y(float y) override { m_state.set_y(y); }

  float get_z() const override { return m_state.get_z(); }
  void set_z(float z) override { m_state.set_z(z); }

  void set_phi(const float phi) override { m_state.set_phi(phi); }
  void set_theta(const float theta) override { m_state.set_theta(theta); }
  void set_qOp(const float qop) override { m_state.set_qOp(qop); }

  float get_pos(unsigned int i) const override { return m_state.get_pos(i); }

  float get_px() const override { return m_state.get_px(); }
  float get_py() const override { return m_state.get_py(); }
  float get_pz() const override { return m_state.get_pz(); }
  float get_mom(unsigned int i) const override { return m_state.get_mom(i); }

  float get_p() const override { return m_state.get_p(); }
  float get_pt() const override { return m_state.get_pt(); }
  float get_eta() const override { return m_state.get_eta(); }
  float get_phi() const override { return m_state.get_phi(); }
  float get_theta() const override { return m_state.get_theta(); }
  float get_qOp() const override { return m_state.get_qOp(); }
  int get_charge() const override { return m_state.get_charge(); }

  float get_covariance(int i, int j) const override { return m_state.get_covariance(i, j); }
  void set_covariance(int i, int j, float value) override { m_state.set_covariance(i, j, value); }

 private:
  // track information
  unsigned int _track_id = std::numeric_limits<unsigned int>::quiet_NaN();
  // unsigned int _vertex_id = UINT_MAX;
  float m_chisq = std::numeric_limits<float>::quiet_NaN();
  uint8_t m_ndf = std::numeric_limits<uint8_t>::quiet_NaN();
  uint64_t m_hitbitmap = std::numeric_limits<uint64_t>::quiet_NaN();
  short int m_crossing = std::numeric_limits<short int>::quiet_NaN();

  // track state information
  TrackStateInfo_v1 m_state;  //< path length => state object

  ClassDefOverride(SvtxTrackInfo_v1, 1)
};

#endif