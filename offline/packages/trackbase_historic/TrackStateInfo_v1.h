#ifndef TRACKBASEHISTORIC_TRACKSTATEINFOV1_H
#define TRACKBASEHISTORIC_TRACKSTATEINFOV1_H

#include <trackbase/TrkrDefs.h>
#include "TrackStateInfo.h"

#include <cmath>
#include <cstddef>  // for size_t
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>  // for pair

class PHObject;

class TrackStateInfo_v1 : public TrackStateInfo
{
 public:
  TrackStateInfo_v1() = default;

  //* base class copy constructor
  // TrackStateInfo_v1( const TrackStateInfo& ) {}

  //* copy constructor
  // TrackStateInfo_v1(const TrackStateInfo_v1& ) {}

  //* assignment operator
  // TrackStateInfo_v1& operator=(const TrackStateInfo_v1& source);

  //* destructor
  ~TrackStateInfo_v1() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackStateInfo_v1 class" << std::endl;
  }
  void Reset() override { *this = TrackStateInfo_v1(); }
  // int isValid() const override;
  PHObject* CloneMe() const override { return new TrackStateInfo_v1(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content from base class
  void CopyFrom(const TrackStateInfo&) override {}
  void CopyFrom(TrackStateInfo* source) override
  {
    CopyFrom(*source);
  }

  //
  // basic track information ---------------------------------------------------
  //

  float get_x() const override { return m_Position[0]; }
  void set_x(float x) override { m_Position[0] = x; }

  float get_y() const override { return m_Position[1]; }
  void set_y(float y) override { m_Position[1] = y; }

  float get_z() const override { return m_Position[2]; }
  void set_z(float z) override { m_Position[2] = z; }

  float get_pos(unsigned int i) const override { return m_Position[i]; }

  float get_px() const override { return m_Momentum[0]; }
  void set_px(float px) override { m_Momentum[0] = px; }

  float get_py() const override { return m_Momentum[1]; }
  void set_py(float py) override { m_Momentum[1] = py; }

  float get_pz() const override { return m_Momentum[2]; }
  void set_pz(float pz) override { m_Momentum[2] = pz; }

  float get_mom(unsigned int i) const override { return m_Momentum[i]; }

  float get_p() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const override { return asinh(get_pz() / get_pt()); }
  float get_phi() const override { return atan2(get_py(), get_px()); }

  // float get_error(int i, int j) const override { return _states.find(0.0)->second->get_error(i, j); }
  // void set_error(int i, int j, float value) override { return _states[0.0]->set_error(i, j, value); }

  float get_covariance(int i, int j) const override;
  void set_covariance(int i, int j, float value) override;

 private:
  float m_Momentum[3] = {std::numeric_limits<float>::quiet_NaN()};  //[-100,100,16] //[x,y,z]
  float m_Position[3] = {std::numeric_limits<float>::quiet_NaN()};  //[-30,30,20]  //[px,py,pz]
  float m_Covariance[21] = {std::numeric_limits<float>::quiet_NaN()};

  // Use m_Covariance[21] instead of [15] for now
  // later on may convert to rotated and then cut to 5X5

  ClassDefOverride(TrackStateInfo_v1, 1)
};

#endif