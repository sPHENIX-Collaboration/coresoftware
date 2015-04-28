#ifndef __EVENT_IO_H__
#define __EVENT_IO_H__

#ifndef LVL2_WINNT
#include "ioselect.h"
#else
#define STREAMBUF_NEW_IOSTREAM
#endif

#if defined(LVL2_WINNT) || defined(STREAMBUF_NEW_IOSTREAM)
#include <iostream>
#include <streambuf>
#include <iomanip>

#ifndef __CINT__

#define COUT std::cout
#define OSTREAM std::ostream
#define ISTREAM std::istream
#define SETW std::setw
#define SETFILL std::setfill
#define ENDL std::endl
#define STREAMBUF std::streambuf

#endif

#include <stdio.h>

#else     // not  LVL2_WINNT

#include <iostream.h>
#include <iomanip.h>

#ifndef __CINT__

#define COUT cout
#define OSTREAM ostream
#define ISTREAM istream
#define SETW setw
#define SETFILL setfill
#define ENDL endl
#define STREAMBUF streambuf

#endif  

#include <stdio.h>

#endif


#endif
