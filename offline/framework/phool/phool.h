#ifndef __PHOOL_H__
#define __PHOOL_H__

//  Standard PHOOL's header file.
//  Purpose: declarations and definitions for PHOOL

#include <iostream>

static const int False = 0;
static const int True = 1;

// PHENIX
static const unsigned int EAST = 0;  // negative x
static const unsigned int WEST = 1;  // positive x
static const unsigned int SOUTH = 0;  // negative z
static const unsigned int NORTH = 1;  // positive z

//  Global type definitions
typedef int PHBoolean;
enum PHMessageType {PHError, PHWarning, PHHullo};
enum PHAccessType {PHReadOnly, PHWrite, PHUpdate};
enum PHTreeType {PHEventTree, PHRunTree};

//  Constants
const double Pi           = 3.14159265358979323846264338328;
const double TwoPi        = 2.0 * Pi;
const double ToRadian     = Pi / 180.0;
const double ToDegree     = 180.0 / Pi;

// General purpose functions
class PHString;
void PHMessage(const PHString&, int, const PHString&);

#ifndef __CINT__
#define PHWHERE __FILE__ << ":" << __LINE__ << ": "
#define PHMESSAGE(x) do {std::cout << PHWHERE << (x) << std::endl;} while(0)
#define PHOOL_VIRTUAL_WARNING do {std::cout << PHWHERE << "using virtual function, doing nothing" << std::endl;} while (0)
// now one where you can give an argument, e.g. the method name
#define PHOOL_VIRTUAL_WARN(x) do {std::cout << PHWHERE << "using virtual function " << x << " doing nothing" << std::endl;} while (0)
#endif



#endif /* __PHOOL_H__ */
