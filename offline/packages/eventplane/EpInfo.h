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
  virtual void InitializeToZero();
  virtual void AddToQraw(unsigned int order, unsigned int index, double val) {QrawOneSide[order][index] += val;}
  virtual void SetWheelSumWeightsRaw(unsigned int order, double val) {WheelSumWeightsRaw[order] = val;}
  virtual void CopyQrawOneSide(const std::vector<std::vector<double>> &vecvec) {QrawOneSide = vecvec;}

 private:
  bool ArgumentOutOfBounds(unsigned int order);

  double Range(double psi, unsigned int order);  /// puts angle psi into range (0,2pi/n)

  std::vector<std::vector<double>> QrawOneSide;
  double WheelSumWeightsRaw[_EpOrderMax];  /// indices: [order]

  ClassDefOverride(EpInfo, 1);
};

#endif  // EVENTPLANE_EPINFO_H
