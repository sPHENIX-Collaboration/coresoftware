#include "SvtxAlignmentState.h"

namespace
{
  SvtxAlignmentState::GlobalMatrix globalMatrix = SvtxAlignmentState::GlobalMatrix::Zero();
  SvtxAlignmentState::LocalMatrix localMatrix = SvtxAlignmentState::LocalMatrix::Zero();
  SvtxAlignmentState::LocalMatrixPsuedo localMatrixPsuedo = SvtxAlignmentState::LocalMatrixPsuedo::Zero();
  SvtxAlignmentState::LocalMeasPsuedoEigenVal localMeasPsuedoEigenVal = SvtxAlignmentState::LocalMeasPsuedoEigenVal::Zero();
  SvtxAlignmentState::ResidualVector residual = SvtxAlignmentState::ResidualVector::Zero();
  SvtxAlignmentState::ActsTrackParamsVector trackParams = SvtxAlignmentState::ActsTrackParamsVector::Zero();
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

const SvtxAlignmentState::LocalMeasPsuedoEigenVal& SvtxAlignmentState::get_local_eigenvalue_psuedo_measurement() const
{
  return localMeasPsuedoEigenVal;
}

const SvtxAlignmentState::GlobalMatrix& SvtxAlignmentState::get_global_derivative_matrix() const
{
  return globalMatrix;
}

const SvtxAlignmentState::ActsTrackParamsVector& SvtxAlignmentState::get_acts_track_params() const
{
  return trackParams;
}
