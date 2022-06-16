#ifndef EVENTPLANE_EPINFO_H
#define EVENTPLANE_EPINFO_H

#define _EpOrderMax 3

#include <phool/PHObject.h>

#include <map>

class EpInfo : public PHObject
{
  friend class EpFinder;

 public:
  EpInfo();
  ~EpInfo() override {/* no op */};

  void Reset() override;

  virtual double RawPsi(unsigned int order);

  virtual double SWRaw(unsigned int order);

  double PsiRaw[_EpOrderMax];

  virtual std::pair<double, double> RawQ(unsigned int order);

  virtual unsigned int MaxOrder() const {return  _EpOrderMax;}

 private:
  bool ArgumentOutOfBounds(unsigned int order);

  double Range(double psi, unsigned int order);  /// puts angle psi into range (0,2pi/n)

  double QrawOneSide[_EpOrderMax][2];      /// indices: [order][x,y]
  double WheelSumWeightsRaw[_EpOrderMax];  /// indices: [order]

  ClassDefOverride(EpInfo, 1);
};

#endif  // EVENTPLANE_EPINFO_H
