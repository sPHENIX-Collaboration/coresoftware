// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHBBOX_H
#define G4MAIN_PHBBOX_H

#include <iostream>
#include <iomanip>
#include <bitset>

// Clip a line using the Cohen-Southerland algorithm

class PHBBox 
{
public:

  //! Construct with the boundaries in x and y
  PHBBox(const double x0, const double y0, const double x1, const double y1) :
    _x0(x0), _y0(y0), _x1(x1), _y1(y1) 
  {}

  //! Given a line, clip it.  Return true if the line should appear in the
  //! bounding box.  The endpoints are updated to reflect the clipping
  bool ClipLine(double& x0, double& y0, double& x1, double& y1) const
  {
    int clipCode0 = ClipCode(x0,y0); // clipping code for end point 0
    int clipCode1 = ClipCode(x1,y1); // clipping code for end point 1

    while ( clipCode0 || clipCode1 )
      {
	//std::cout << "clipCode0 = " << std::bitset<4>(clipCode0).to_string() << std::endl;
	//std::cout << "clipCode1 = " << std::bitset<4>(clipCode1).to_string() << std::endl;

	if ( clipCode0 & clipCode1 ) return false;

	int code = 0;
	if ( clipCode0 > 0 ) code = clipCode0; // clip the first point
	else                 code = clipCode1; // clip the last point

	double x = 0, y = 0;

	if ( (code & BOTTOM) == BOTTOM )
	  {
	    // Clip the line to the bottom of the box
	    //std::cout << "Clip the line to the bottom of the box" << std::endl;
	    y = _y0;
	    x = x0 + (x1-x0)*(y-y0)/(y1-y0);
	  }
	else if ( (code & TOP) == TOP ) 
	  {
	    // Clip the line to the top of the box
	    //std::cout << "Clip the line to the top of the box" << std::endl;
	    y = _y1;
	    x = x0 + (x1-x0)*(y-y0)/(y1-y0);
	  }
	else if ( (code & LEFT) == LEFT )
	  {
	    //std::cout << "Clip the line to the left of the box" << std::endl;
	    x = _x0;
	    y = y0 + (y1-y0)*(x-x0)/(x1-x0);
	  }
	else if ( (code & RIGHT) == RIGHT )
	  {
	    //std::cout << "Clip the line to the right of the box" << std::endl;
	    x = _x1;
	    y = y0 + (y1-y0)*(x-x0)/(x1-x0);
	  }

	//std::cout << "x = " << x << ", y = " << y << std::endl;

	if ( code == clipCode0 )
	  {
	    // modify the first coord
	    //std::cout << "modify the first coord" << std::endl;
	    x0 = x;
	    y0 = y;
	    clipCode0 = ClipCode(x0,y0);
	  }
	else 
	  {
	    // modify the second coord
	    //std::cout << "modify the second coord" << std::endl;
	    x1 = x;
	    y1 = y;
	    clipCode1 = ClipCode(x1,y1);
	  }
      }

    return true;
  }

  void Print(std::ostream& os = std::cout)
  {
    os << _x0 << " " << _y0 << ", " << _x1 << " " << _y1 << std::endl;
  }

private:
  enum { RIGHT=1, BOTTOM=2, LEFT=4, TOP=8 };

  int ClipCode(const double x, const double y) const
  {
    int code = 0;
    if      ( x > _x1 ) code |= RIGHT;
    else if ( x < _x0 ) code |= LEFT;
    if      ( y > _y1 ) code |= TOP;
    else if ( y < _y0 ) code |= BOTTOM;
    return code;
  }

  double _x0;
  double _y0;
  double _x1;
  double _y1; 
};

#endif // __PHBBOX_H__
