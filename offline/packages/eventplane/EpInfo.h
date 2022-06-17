#ifndef EVENTPLANE_EPINFO_H
#define EVENTPLANE_EPINFO_H

#include <phool/PHObject.h>

#include <cmath>
#include <map>

class EpInfo : public PHObject
{
 public:
  EpInfo() {}
  ~EpInfo() override {}

  void Reset() override { return; }

  virtual double RawPsi(unsigned int /* order */) const { return NAN; }

  virtual double SWRaw(unsigned int /* order */) const { return NAN; }

  virtual std::pair<double, double> RawQ(unsigned int /* order */) const { return std::make_pair(NAN, NAN); }

  virtual unsigned int MaxOrder() const { return 0; }
  virtual void CopyQrawOneSide(const std::vector<std::vector<double>>& /*vecvec */) { return; }
  virtual void CopyWheelSumWeightsRaw(const std::vector<double>& /*vec */) { return; }
  virtual void CopyPsiRaw(const std::vector<double>& /*vec */) { return; }

 private:
  ClassDefOverride(EpInfo, 1);
};

#endif  // EVENTPLANE_EPINFO_H
