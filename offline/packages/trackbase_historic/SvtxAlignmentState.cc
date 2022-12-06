#include "SvtxAlignmentState.h"

namespace
{
  SvtxAlignmentState::GlobalMatrix globalMatrix = SvtxAlignmentState::GlobalMatrix::Zero();
  SvtxAlignmentState::LocalMatrix localMatrix = SvtxAlignmentState::LocalMatrix::Zero();
  SvtxAlignmentState::ResidualVector residual = SvtxAlignmentState::ResidualVector::Zero();
}  // namespace

const SvtxAlignmentState::ResidualVector& SvtxAlignmentState::get_residual() const
{
  return residual;
}

const SvtxAlignmentState::LocalMatrix& SvtxAlignmentState::get_local_derivative_matrix() const
{
  return localMatrix;
}

const SvtxAlignmentState::GlobalMatrix& SvtxAlignmentState::get_global_derivative_matrix() const
{
  return globalMatrix;
}
