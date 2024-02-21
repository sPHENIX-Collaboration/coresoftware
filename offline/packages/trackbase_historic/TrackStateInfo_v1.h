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

  ~TrackStateInfo_v1() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackStateInfo_v1 class" << std::endl;
  }
  void Reset() override { *this = TrackStateInfo_v1(); }

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

  float get_px() const override;
  float get_py() const override;
  float get_pz() const override;

  float get_phi() const override { return m_Momentum[0]; }
  void set_phi(const float phi) override { m_Momentum[0] = phi; }

  int get_charge() const override { return get_qOp() > 0 ? 1 : -1; }
  float get_theta() const override { return m_Momentum[1]; }
  void set_theta(const float theta) override { m_Momentum[1] = theta; }

  float get_qOp() const override { return m_Momentum[2]; }
  void set_qOp(const float qop) override { m_Momentum[2] = qop; }

  float get_mom(unsigned int i) const override
  {
    if (i == 0)
    {
      return get_px();
    }
    if (i == 1)
    {
      return get_py();
    }
    if (i == 2)
    {
      return get_pz();
    }
    return std::numeric_limits<float>::quiet_NaN();
  }

  float get_p() const override { return sqrt(std::pow(get_px(), 2) + std::pow(get_py(), 2) + std::pow(get_pz(), 2)); }
  float get_pt() const override { return sqrt(std::pow(get_px(), 2) + std::pow(get_py(), 2)); }
  float get_eta() const override { return -std::log(std::tan(get_theta() / 2.)); }

  float get_covariance(int i, int j) const override;
  void set_covariance(int i, int j, float value) override;

 private:
  float m_Momentum[3] = {std::numeric_limits<float>::quiet_NaN()};  // global phi, theta, q/p
  float m_Position[3] = {std::numeric_limits<float>::quiet_NaN()};  // global [x,y,z]
  float m_Covariance[15] = {std::numeric_limits<float>::quiet_NaN()};

  ClassDefOverride(TrackStateInfo_v1, 1)
};

#endif