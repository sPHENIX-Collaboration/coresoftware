// @sa <O2/DataFormats/common/src/InteractionRecord.cxx>
//     <03608ff89>

#include "mvtx_decoder/InteractionRecord.h"
#include <iostream>
#include <sstream>

namespace mvtx
{

  //________________________________________________________________________________
  std::string InteractionRecord::asString() const
  {
    if (isDummy())
    {
      return std::string{"NotSet"};
    }

    std::stringstream ss;
    ss << std::hex << "Orbit: 0x" << orbit << std::dec << "  BCid: " << bc;
    return ss.str();
  }

  //________________________________________________________________________________
  std::ostream& operator<<(std::ostream& stream, mvtx::InteractionRecord const& ir)
  {
    stream << ir.asString();
    return stream;
  }

  //________________________________________________________________________________
  void InteractionRecord::print() const
  {
    std::cout << (*this) << std::endl;
  }

}  // namespace mvtx
