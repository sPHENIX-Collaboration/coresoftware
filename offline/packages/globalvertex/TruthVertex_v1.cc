#include "TruthVertex_v1.h"

void TruthVertex_v1::identify(std::ostream& os) const
{
  os << "TruthVertex_v1 - ID: " << _id
     << ", z: " << _z << " ± " << _z_err
     << ", t: " << _t << " ± " << _t_err
     << std::endl;
}

int TruthVertex_v1::isValid() const
{
  return std::isfinite(_z) && std::isfinite(_t);
}
