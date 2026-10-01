#include "SvtxAlignmentState.h"

namespace
{
  SvtxAlignmentState::GlobalMatrix globalMatrix = SvtxAlignmentState::GlobalMatrix::Zero();
  SvtxAlignmentState::LocalMatrix localMatrix = SvtxAlignmentState::LocalMatrix::Zero();
  SvtxAlignmentState::LocalMatrixPsuedo localMatrixPsuedo = SvtxAlignmentState::LocalMatrixPsuedo::Zero();
  SvtxAlignmentState::ActsTrackParamsVector localMeasErrPsuedo = SvtxAlignmentState::ActsTrackParamsVector::Zero();
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

const SvtxAlignmentState::LocalMatrixPsuedo& SvtxAlignmentState::get_local_derivative_psuedo_matrix() const
{
  return localMatrixPsuedo;
}

const SvtxAlignmentState::ActsTrackParamsVector& SvtxAlignmentState::get_local_psuedo_measurement_err() const
{
  return localMeasErrPsuedo;
}

const SvtxAlignmentState::GlobalMatrix& SvtxAlignmentState::get_global_derivative_matrix() const
{
  return globalMatrix;
}

