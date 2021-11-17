#include "EventHeader.h"
#include <iostream>

/// Clear Event
void EventHeader::Reset()
{
  std::cout << __PRETTY_FUNCTION__
            << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

/** identify Function from PHObject
 @param os Output Stream
 */
void EventHeader::identify(std::ostream& os) const
{
  os << "identify yourself: virtual EventHeader Object" << std::endl;
  return;
}

/// isValid returns non zero if object contains valid data
int EventHeader::isValid() const
{
  std::cout << __PRETTY_FUNCTION__ << "isValid not implemented by daughter class"
            << std::endl;
  return 0;
}
