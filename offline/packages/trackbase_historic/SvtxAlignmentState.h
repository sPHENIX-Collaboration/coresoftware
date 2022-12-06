#ifndef __SVTXALIGNMENTSTATE_H__
#define __SVTXALIGNMENTSTATE_H__

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <Acts/Definitions/Algebra.hpp>
#include <cmath>

class SvtxAlignmentState : public PHObject
{
 public:
  /// Number of global coordinates, where global refers to
  /// the definition in millepede (alignment parameters)
  const static int NGL = 6;
  /// Number of local coordinates, where local refers to
  /// the definition in millepede (track state parameters)
  const static int NLOC = 6;
  /// Number of residual parameters
  const static int NRES = 3;

  typedef Eigen::Matrix<double, NRES, NGL> GlobalMatrix;
  typedef Eigen::Matrix<double, NRES, NLOC> LocalMatrix;
  typedef Eigen::Matrix<double, NRES, 1> ResidualVector;

  ~SvtxAlignmentState() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxAlignmentState base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual void set_residual(const ResidualVector&) {}
  virtual void set_local_derivative_matrix(const LocalMatrix&) {}
  virtual void set_global_derivative_matrix(const GlobalMatrix&) {}
  virtual void set_cluster_key(TrkrDefs::cluskey) {}

  virtual const ResidualVector& get_residual() const;
  virtual const LocalMatrix& get_local_derivative_matrix() const;
  virtual const GlobalMatrix& get_global_derivative_matrix() const;
  virtual TrkrDefs::cluskey get_cluster_key() const { return UINT_MAX; }

 protected:
  SvtxAlignmentState() {}
  ClassDefOverride(SvtxAlignmentState, 1);
};

#endif
