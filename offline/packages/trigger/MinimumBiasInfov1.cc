#include "MinimumBiasInfov1.h"

#include <ostream>

void MinimumBiasInfov1::identify(std::ostream& os) const
{
<<<<<<< HEAD

  os << "MinimumBiasInfo: "<< std::endl;
  os << "  IsMinBias = " << (_isMinBias? "Yes":"No") << std::endl;
=======
  os << "MinimumBiasInfo: " << std::endl;
  os << "  IsMinBias = " << (_isMinBias ? "Yes" : "No") << std::endl;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return;
}
