#ifndef RHOUECALIB_MBDCHARGEHISTMAKER_H
#define RHOUECALIB_MBDCHARGEHISTMAKER_H

#include "RhoUECalibBase.h"

#include <string>

class PHCompositeNode;

/// \class MbdChargeHistMaker
///
/// \brief Step 1 of the rho UE calibration: the MBD charge Q = q_S + q_N of the selected
/// events (minimum bias, vertex in the class z range), one histogram per run
/// ("q_<run>", 1 unit bins up to set_q_max), and the MBD t0 per run ("t0_<run>", 0.05 ns
/// bins, |t0| < 25 ns). RhoUECalibPayload::FitMbdKnees fits the central edge ("knee") of each
/// charge spectrum and finds the main-bunch centre of each t0 spectrum, and writes the knee
/// table of the later steps (RhoUECalibBase::set_knee_table, RhoUECalibEventVeto::set_t0_table).
class MbdChargeHistMaker : public RhoUECalibBase
{
 public:
  explicit MbdChargeHistMaker(const std::string &outfile, const std::string &name = "MbdChargeHistMaker")
    : RhoUECalibBase(name, outfile)
  {
  }

  void set_q_max(const double qmax) { m_qmax = qmax; }

 protected:
  int fill(PHCompositeNode *topNode) override;
  void book() override {}
  bool needs_class() const override { return false; }
  bool needs_geometry() const override { return false; }

 private:
  double m_qmax{4000};
};

#endif
