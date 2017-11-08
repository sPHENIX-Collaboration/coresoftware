#ifndef __PHOOL_H__
#define __PHOOL_H__

//  Standard PHOOL's header file.
//  Purpose: declarations and definitions for PHOOL

#include <iostream>

static const int False = 0;
static const int True = 1;

//  Global type definitions
typedef int PHBoolean;
enum PHMessageType {PHError, PHWarning, PHHullo};
enum PHAccessType {

  //! Read from DST file
  PHReadOnly,

  //! Write to DST file
  PHWrite,

  //! Update DST file
  PHUpdate,

  //! Read from DST file and integrate its content with the current nodes
  PHReadAndIntegrate
};
enum PHTreeType {
  //! DST node, which store event-wise information
  PHEventTree = 0,

  //! RUN node, which store run-wise information
  PHRunTree = 1,

  //! SUM node,  which store integrals, e.g. integrated luminosity
  PHIntegralTree = 2
};

// General purpose functions
void PHMessage(const std::string&, int, const std::string&);

#ifndef __CINT__
#define PHWHERE __FILE__ << ":" << __LINE__ << ": "
#define PHMESSAGE(x) do {std::cout << PHWHERE << (x) << std::endl;} while(0)
#define PHOOL_VIRTUAL_WARNING do {std::cout << PHWHERE << "using virtual function, doing nothing" << std::endl;} while (0)
// now one where you can give an argument, e.g. the method name
#define PHOOL_VIRTUAL_WARN(x) do {std::cout << PHWHERE << "using virtual function " << x << " doing nothing" << std::endl;} while (0)
#endif



#endif /* __PHOOL_H__ */
