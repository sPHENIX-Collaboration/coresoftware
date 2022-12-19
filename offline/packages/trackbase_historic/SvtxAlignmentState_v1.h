#ifndef TRACKBASEHISTORIC_SVTXALIGNMENTSTATE_V1_H
#define TRACKBASEHISTORIC_SVTXALIGNMENTSTATE_V1_H

#include "SvtxAlignmentState.h"

#include <cmath>
#include <iostream>

class PHObject;

class SvtxAlignmentState_v1 : public SvtxAlignmentState
{
 public:
  SvtxAlignmentState_v1();
  ~SvtxAlignmentState_v1() override {}

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxAlignmentState_v1(); }
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxAlignmentState_v1(*this); }

  void set_residual(const ResidualVector& res) override
  {
    m_residual = res;
  }
  void set_local_derivative_matrix(const LocalMatrix& d) override
  {
    m_localDeriv = d;
  }
  void set_global_derivative_matrix(const GlobalMatrix& d) override
  {
    m_globalDeriv = d;
  }
  void set_cluster_key(const TrkrDefs::cluskey key) override
  {
    m_cluskey = key;
  }

  const ResidualVector& get_residual() const override { return m_residual; }
  const LocalMatrix& get_local_derivative_matrix() const override { return m_localDeriv; }
  const GlobalMatrix& get_global_derivative_matrix() const override { return m_globalDeriv; }
  TrkrDefs::cluskey get_cluster_key() const override { return m_cluskey; }

 private:
  ResidualVector m_residual;
  LocalMatrix m_localDeriv;
  GlobalMatrix m_globalDeriv;
  TrkrDefs::cluskey m_cluskey;

  ClassDefOverride(SvtxAlignmentState_v1, 1)
};

#endif
