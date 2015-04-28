#ifndef __BOOLDEFS_H
#define __BOOLDEFS_H

/*
 *--------------------------------------------
 * 
 * booldefs.h
 *
 *--------------------------------------------
 */


/**
 Problem: how to handle bool data type when the
 CC compiler doesn't know it?.
 Tricky example:
 Irix6: on R10000 CC has it built in, but if
 compiled with -o32 it is not there.
 In addition STL defines a bool that will give
 a conflict if you define it by hand.


 @author M.Purschke, C.Witzig \\
 {\bf Date: } April 8, 98 \\

 @version Last update April 8, 98

*/


#ifndef BOOL
#define BOOL int
#define true 1
#define false  0
#endif


#endif



