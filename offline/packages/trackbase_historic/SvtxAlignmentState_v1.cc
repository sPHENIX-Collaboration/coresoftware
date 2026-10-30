#include "SvtxAlignmentState_v1.h"

SvtxAlignmentState_v1::SvtxAlignmentState_v1()
  : m_residual(ResidualVector::Zero())
  , m_localDeriv(LocalMatrix::Zero())
  , m_globalDeriv(GlobalMatrix::Zero())
  , m_localDerivPsuedo(LocalMatrixPsuedo::Zero())
  , m_localMeasErrPsuedo(ActsTrackParamsVector::Zero())
  , m_cluskey(UINT_MAX)
{
}

void SvtxAlignmentState_v1::identify(std::ostream &os) const
{
  os << "SvtxAlignmentState_v1 identify: " << std::endl;
  os << "state for cluskey : " << m_cluskey << std::endl;
  os << "residual : " << m_residual.transpose() << std::endl;
  os << "local derivatives : " << std::endl
     << m_localDeriv << std::endl;
  os << "global derivatives : " << std::endl
     << m_globalDeriv << std::endl;
  os << "local pseudo-measurement derivatives : " << std::endl
     << m_localDerivPsuedo << std::endl;
  os << "local pseudo-measurement errors : " << m_localMeasErrPsuedo.transpose() << std::endl;
}
