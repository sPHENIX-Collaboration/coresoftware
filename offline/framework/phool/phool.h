#ifndef PHOOL_PHOOL_H
#define PHOOL_PHOOL_H

//  Standard PHOOL's header file.
//  Purpose: declarations and definitions for PHOOL

#include <iostream>

//  Global type definitions
enum PHMessageType
{
  PHError,
  PHWarning,
  PHHullo
};
enum PHAccessType
{
  PHReadOnly,
  PHWrite,
  PHUpdate
};
enum PHTreeType
{
  PHEventTree,
  PHRunTree
};

// General purpose functions
void PHMessage(const std::string&, int, const std::string&);

#define PHWHERE __FILE__ << ":" << __LINE__ << ": "
#define PHMESSAGE(x)                          \
  do                                          \
  {                                           \
    std::cout << PHWHERE << (x) << std::endl; \
  } while (0)
#define PHOOL_VIRTUAL_WARNING                                                     \
  do                                                                              \
  {                                                                               \
    std::cout << PHWHERE << "using virtual function, doing nothing" << std::endl; \
  } while (0)
// now one where you can give an argument, e.g. the method name
#define PHOOL_VIRTUAL_WARN(x)                                                                \
  do                                                                                         \
  {                                                                                          \
    std::cout << PHWHERE << "using virtual function " << x << " doing nothing" << std::endl; \
  } while (0)

#endif /* PHOOL_PHOOL_H */
