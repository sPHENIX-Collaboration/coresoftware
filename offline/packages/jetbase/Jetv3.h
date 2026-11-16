/*!
 * \file    Jetv3.h
 * \brief   Jetv2 plus the corrected jet pT (MC JES calibration, in-situ JES on data, in-situ
 *          JER smearing on simulation), so that the correction modules write into the jet
 *          itself instead of creating a new jet node.
 */

#ifndef JETBASE_JETV3_H
#define JETBASE_JETV3_H

#include "Jetv2.h"

#include <iostream>
#include <limits>

class PHObject;

/*!
 * \brief Jetv3
 *
 * get_pt() stays the reconstructed (raw) pT. The corrected values are separate members:
 *   pt_calib               after the MC JES calibration (JetCalib)
 *   pt_insitu              after the in-situ JES calibration (data)
 *   pt_insitu_up / _down   its systematic variations
 *   pt_smear               after the in-situ JER smearing (simulation)
 *   pt_smear_up / _down    its systematic variations
 * All are NaN until a module fills them. Stored in JetContainerv2.
 */
class Jetv3 : public Jetv2
{
 public:
  Jetv3() = default;
  explicit Jetv3(unsigned int n_prop);

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  PHObject* CloneMe() const override;

  // corrected jet pT
  float get_pt_calib() const override { return _pt_calib; }
  void set_pt_calib(float pt) override { _pt_calib = pt; }

  float get_pt_insitu() const override { return _pt_insitu; }
  void set_pt_insitu(float pt) override { _pt_insitu = pt; }

  float get_pt_insitu_up() const override { return _pt_insitu_up; }
  void set_pt_insitu_up(float pt) override { _pt_insitu_up = pt; }

  float get_pt_insitu_down() const override { return _pt_insitu_down; }
  void set_pt_insitu_down(float pt) override { _pt_insitu_down = pt; }

  float get_pt_smear() const override { return _pt_smear; }
  void set_pt_smear(float pt) override { _pt_smear = pt; }

  float get_pt_smear_up() const override { return _pt_smear_up; }
  void set_pt_smear_up(float pt) override { _pt_smear_up = pt; }

  float get_pt_smear_down() const override { return _pt_smear_down; }
  void set_pt_smear_down(float pt) override { _pt_smear_down = pt; }

 private:
  float _pt_calib{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu_up{std::numeric_limits<float>::quiet_NaN()};
  float _pt_insitu_down{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear_up{std::numeric_limits<float>::quiet_NaN()};
  float _pt_smear_down{std::numeric_limits<float>::quiet_NaN()};

  ClassDefOverride(Jetv3, 1);
};

#endif  // JETBASE_JETV3_H
