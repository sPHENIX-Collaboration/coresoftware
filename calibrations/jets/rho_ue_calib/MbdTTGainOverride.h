#ifndef RHOUECALIB_MBDTTGAINOVERRIDE_H
#define RHOUECALIB_MBDTTGAINOVERRIDE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

/// \class MbdTTGainOverride
///
/// \brief Sets the MBD time-channel gain (TTGain) of every PMT to one value after
/// MbdReco::InitRun loaded the calibrations (register it directly after MbdReco). Without an
/// MBD_TTGAIN payload in the CDB, MbdCalib falls back to 0.95 for every PMT; the Run 2024
/// Au+Au production timing (t0, z) is reproduced with 1.0. MbdEvent reads the gain per event.
class MbdTTGainOverride : public SubsysReco
{
 public:
  explicit MbdTTGainOverride(const float gain = 1.0, const std::string &mbdreco = "MbdReco",
                             const std::string &name = "MbdTTGainOverride")
    : SubsysReco(name)
    , m_mbdreco(mbdreco)
    , m_gain(gain)
  {
  }
  ~MbdTTGainOverride() override = default;
  MbdTTGainOverride(const MbdTTGainOverride &) = delete;
  MbdTTGainOverride &operator=(const MbdTTGainOverride &) = delete;
  MbdTTGainOverride(MbdTTGainOverride &&) = delete;
  MbdTTGainOverride &operator=(MbdTTGainOverride &&) = delete;

  int InitRun(PHCompositeNode *topNode) override;

 private:
  std::string m_mbdreco;
  float m_gain;
};

#endif
