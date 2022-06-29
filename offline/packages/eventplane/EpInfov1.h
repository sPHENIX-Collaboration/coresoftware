#ifndef EVENTPLANE_EPINFOV1_H
#define EVENTPLANE_EPINFOV1_H

#include "EpInfo.h"

#include <utility>  // for pair
#include <vector>   // for vector

class EpInfov1 : public EpInfo
{
 public:
  EpInfov1() {}
  ~EpInfov1() override{/* no op */};

  void Reset() override;

  virtual double RawPsi(unsigned int order) const override;

  virtual double SWRaw(unsigned int order) const override;

  virtual std::pair<double, double> RawQ(unsigned int order) const override;

  virtual unsigned int MaxOrder() const override { return QrawOneSide.size(); }
  virtual void CopyQrawOneSide(const std::vector<std::vector<double>> &vecvec) override;  // {QrawOneSide = vecvec;}
  virtual void CopyWheelSumWeightsRaw(const std::vector<double> &vec) override;           // {WheelSumWeightsRaw = vec;}
  virtual void CopyPsiRaw(const std::vector<double> &vec) override;                       // {PsiRaw = vec;}

 private:
  bool ArgumentOutOfBounds(unsigned int order) const;

  double Range(double psi, unsigned int order) const;  /// puts angle psi into range (0,2pi/n)

  std::vector<std::vector<double>> QrawOneSide;
  std::vector<double> WheelSumWeightsRaw;
  std::vector<double> PsiRaw;

  ClassDefOverride(EpInfov1, 1);
};

#endif  // EVENTPLANE_EPINFOV1_H
