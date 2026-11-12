#ifndef JETBACKGROUND_MBDCHARGEKNEE_H
#define JETBACKGROUND_MBDCHARGEKNEE_H

#include <string>

class TH1;

/// \class MbdChargeKnee
///
/// \brief MBD charge knee Q_knee of a run, the scale of Q' = (q_S + q_N) / Q_knee in the rho
/// UE calibrations (DetermineTowerRho tower shape map, TowerRhoToBackground)
///
/// Q_knee is the central edge of the run's minimum-bias MBD charge spectrum
/// (RhoUECalibPayload::FitMbdKnees). It is taken, in order, from
///   - a value set by the user (> 0);
///   - the per-run knee payload, when a CDB tag or a path is set: a CDBTTree with the single
///     float "mbdq_knee" (RhoUECalibPayload::MakeKneePayloads, one payload per run). The path is
///     one payload, or a directory with the payloads of many runs, mbdq_knee_<run>.root;
///   - the per-run table of the calibration payload (TH1 over the run number, the runs the
///     calibration was made from);
///   - the dataset default: the median knee of that table (dataset_default is set).
class MbdChargeKnee
{
 public:
  // Q_knee of the run, -1 if there is none; source says where it came from
  static float find(float user, const std::string &tag, const std::string &path, const TH1 *table, int run,
                    std::string &source, bool &dataset_default);
};

#endif
