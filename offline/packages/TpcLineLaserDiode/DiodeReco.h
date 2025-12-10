#ifndef DIODERECO__H
#define DIODERECO__H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class CDBInterface;
class PHCompositeNode;
class TpcDiodeContainer;

class DiodeReco : public SubsysReco
{
 public:
  DiodeReco(const std::string &name = "DiodeReco");
  ~DiodeReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  double MaxAdc(int n, int low_bin = 0, int high_bin = 9999) const;  // Signal is averaged over "avg" number of bins interatively within the bin range [low_bin,high_bin]
  int MaxBin(int n) const;
  double Integral(int low_bin, int high_bin) const;
  int NAboveThreshold(double upper_thr, double lower_thr) const;
  double PulseWidth(double upper_thr, double lower_thr) const;
  void PedestalCorrected(int low_bin, int high_bin);

 private:
  CDBInterface *m_cdb{nullptr};
  TpcDiodeContainer *diodes{nullptr};
  std::string calibdir;
  std::string m_DiodeContainerName;
  std::vector<double> adc;
};

#endif
