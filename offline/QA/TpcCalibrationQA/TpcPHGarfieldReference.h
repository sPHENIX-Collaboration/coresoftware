#ifndef QA_TPCCALIBRATIONQA_PHGARFIELDREFERENCE_H
#define QA_TPCCALIBRATIONQA_PHGARFIELDREFERENCE_H

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <array>
#include <cmath>
#include <memory>

class PHCompositeNode;
class PHGarfield;

// Initializes the equal-FEE PHGarfield baseline once per run; fills no event histograms.
class TpcPHGarfieldReference : public SubsysReco
{
 public:
  TpcPHGarfieldReference(double field, double pressure, double temperature,
                         double t0, double clock, double step);
  ~TpcPHGarfieldReference() override;

  TpcPHGarfieldReference(const TpcPHGarfieldReference &other) = delete;
  TpcPHGarfieldReference &operator=(const TpcPHGarfieldReference &other) = delete;
  TpcPHGarfieldReference(TpcPHGarfieldReference &&other) = delete;
  TpcPHGarfieldReference &operator=(TpcPHGarfieldReference &&other) = delete;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  const std::array<double, 2> &velocity() const
  {
    return m_velocity;
  }
  double pressure() const
  {
    return m_pressure;
  }
  double temperatureK() const
  {
    return m_temperatureK;
  }
  double effectiveField() const
  {
    return m_effectiveField;
  }
  double padPlane() const
  {
    return m_zPadPlaneCm;
  }
  double cmHalf() const
  {
    return m_cmHalf;
  }

 private:
  std::array<double, 2> m_velocity{};
  double m_cmHalf{0.};
  static constexpr int kNFees = 312;
  std::array<std::array<double, kNFees>, 2> m_feeX{}, m_feeY{}, m_feeCount{};
  std::array<std::array<double, kNFees>, 2> m_phgDriftTimeNs{}, m_phgDriftVelocityCmNs{};
  bool m_usePTCorrectedPHGarfield{true};
  double m_referenceCMFieldVcm, m_referencePressure, m_referenceTemperatureK;
  double m_tpcTimeBinNs, m_phgT0Bins, m_phgStepNs, m_zPadPlaneCm{0.};
  double m_pressure{0.}, m_temperatureK{0.}, m_effectiveField{0.};
  std::unique_ptr<PHGarfield> m_driftPHGarfield;
  static bool is_good_number(double v)
  {
    return std::isfinite(v) && std::fabs(v) < 1.e30;
  }
  static double temperature_kelvin(double v)
  {
    return v < 200. ? v + 273.15 : v;
  }
  bool buildFeeMap(PHCompositeNode *topNode);
  bool initializeDriftCoordinateCalibration(PHCompositeNode *topNode);
};
#endif
