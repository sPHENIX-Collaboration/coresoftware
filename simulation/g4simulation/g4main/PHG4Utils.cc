#include "PHG4Utils.h"

#include <phool/phool.h>

#include <Geant4/G4Colour.hh>  // for G4Colour
#include <Geant4/G4VisAttributes.hh>

#include <TSystem.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/algorithm/string.hpp>
#pragma GCC diagnostic pop

#include <boost/algorithm/hex.hpp>
#include <boost/uuid/detail/md5.hpp>

#include <algorithm>  // for copy
#include <cmath>
#include <cstdlib>  // for exit
#include <fstream>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <iterator>  // for back_insert_iterator
#include <vector>    // for vector

using namespace std;

double PHG4Utils::_eta_coverage = 1.;

double
PHG4Utils::GetLengthForRapidityCoverage(const double radius, const double eta)
{
  double length;
  double theta = 2.0 * std::atan(std::exp(-eta));
  length = radius / std::tan(theta);
  return length;
}

double
PHG4Utils::GetLengthForRapidityCoverage(const double radius)
{
  return GetLengthForRapidityCoverage(radius, _eta_coverage);
}

void PHG4Utils::SetPseudoRapidityCoverage(const double eta)
{
  _eta_coverage = eta;
}

double
PHG4Utils::get_theta(const double eta)
{
  double theta = 2 * atan(exp(-eta));
  return theta;
}

double
PHG4Utils::get_eta(const double theta)
{
  double eta = -log(tan(theta / 2.));
  return eta;
}

pair<double, double>
PHG4Utils::get_etaphi(const double x, const double y, const double z)
{
  double eta;
  double phi;
  double radius;
  double theta;
  radius = sqrt(x * x + y * y);
  phi = atan2(y, x);
  theta = atan2(radius, z);
  eta = -log(tan(theta / 2.));
  return make_pair(eta, phi);
}

double
PHG4Utils::get_eta(const double radius, const double z)
{
  double eta;
  double theta;
  theta = atan2(radius, fabs(z));
  eta = -log(tan(theta / 2.));
  if (z < 0)
  {
    eta = -eta;
  }
  return eta;
}

void PHG4Utils::SetColour(G4VisAttributes* att, const string& material)
{
  if (!att)
  {
    cout << "G4VisAttributes pointer is NULL" << endl;
    return;
  }
  if (material == "AL_BABAR_MAG")
  {
    att->SetColour(G4Colour::Blue());
  }
  else if (material == "BlackHole")
  {
    att->SetColour(G4Colour::Black());
  }
  else if (material == "C4F10")
  {
    att->SetColour(0., 0., 0.5, 0.25);
  }
  else if (material == "CF4")
  {
    att->SetColour(G4Colour::Magenta());
  }
  else if (material == "G4_AIR")
  {
    att->SetColour(G4Colour::Black());
  }
  else if (material == "G4_Al")
  {
    att->SetColour(G4Colour::Gray());
  }
  else if (material == "G4_Au")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_CARBON_DIOXIDE")
  {
    att->SetColour(G4Colour::Green());
  }
  else if (material == "G4_CELLULOSE_CELLOPHANE")
  {
    att->SetColour(0.25, 0.25, 0.);
  }
  else if (material == "G4_Cu")
  {
    att->SetColour(1., 0.51, 0.278);
  }
  else if (material == "G4_Fe")
  {
    att->SetColour(0.29, 0.44, 0.54);
  }
  else if (material == "G4_KAPTON")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_MYLAR")
  {
    att->SetColour(0.5, 0.5, 0.5, 0.25);
  }
  else if (material == "G4_METHANE")
  {
    att->SetColour(0., 1., 1., 0.25);
  }
  else if (material == "G4_Si")
  {
    att->SetColour(G4Colour::Yellow());
  }
  else if (material == "G4_TEFLON")
  {
    att->SetColour(G4Colour::White());
  }
  else if (material == "G4_W")
  {
    att->SetColour(0.36, 0.36, 0.36);
  }
  else if (material == "Quartz")
  {
    att->SetColour(G4Colour::Green());
  }
  else if (material == "Scintillator" || material == "G4_POLYSTYRENE")
  {
    att->SetColour(0., 1., 1.);
  }
  else if (material == "W_Epoxy")
  {
    att->SetColour(0.5, 0.5, 0.5);
  }
  else if (material == "G10")
  {
    att->SetColour(1., 1., 0., 0.5);
  }
  else
  {
    //cout << "default color red for material " << material << endl;
    att->SetColour(G4Colour::Cyan());
  }
  return;
}
// returns success/failure and intersection point (x/y)
std::pair<bool, std::pair<double, double>> PHG4Utils::lines_intersect(
    double ax,
    double ay,
    double bx,
    double by,
    double cx,
    double cy,
    double dx,
    double dy)
{
  // Find if a line segment limited by points A and B
  // intersects line segment limited by points C and D.
  // First check if an infinite line defined by A and B intersects
  // segment (C,D). If h is from 0 to 1 line and line segment intersect
  // Then check in intersection point is between C and D
  double ex = bx - ax;  // E=B-A
  double ey = by - ay;
  double fx = dx - cx;  // F=D-C
  double fy = dy - cy;
  double px = -ey;  // P
  double py = ex;

  double bottom = fx * px + fy * py;  // F*P
  double gx = ax - cx;                // A-C
  double gy = ay - cy;
  double top = gx * px + gy * py;  // G*P

  double h = 99999.;
  if (bottom != 0.)
  {
    h = top / bottom;
  }

  //intersection point R = C + F*h
  if (h > 0. && h < 1.)
  {
    double rx = cx + fx * h;
    double ry = cy + fy * h;
    //cout << "      line/segment intersection coordinates: " << *rx << " " << *ry << endl;
    if ((rx > ax && rx > bx) || (rx < ax && rx < bx) || (ry < ay && ry < by) || (ry > ay && ry > by))
    {
      //cout << "       NO segment/segment intersection!" << endl;
      return make_pair(false, make_pair(NAN, NAN));
    }
    else
    {
      //cout << "       segment/segment intersection!" << endl;
      return make_pair(true, make_pair(rx, ry));
    }
  }

  return make_pair(false, make_pair(NAN, NAN));
}

// returns success/failure and length of the line segment inside the rectangle (output)
std::pair<bool, double> PHG4Utils::line_and_rectangle_intersect(
    double ax,
    double ay,
    double bx,
    double by,
    double cx,
    double cy,
    double dx,
    double dy)
{
  // find if a line isegment limited by points (A,B)
  // intersects with a rectangle defined by two
  // corner points (C,D) two other points are E and F
  //   E--------D
  //   |        |
  //   |        |
  //   C--------F

  if (cx > dx || cy > dy)
  {
    cout << PHWHERE << "ERROR: Bad rectangle definition!" << endl;
    return make_pair(false, NAN);
  }

  double ex = cx;
  double ey = dy;
  double fx = dx;
  double fy = cy;

  vector<double> vx;
  vector<double> vy;
  pair<bool, pair<double, double>> intersect1 = lines_intersect(ax, ay, bx, by, cx, cy, fx, fy);
  //  bool i1 = lines_intersect(ax, ay, bx, by, cx, cy, fx, fy, &rx, &ry);
  if (intersect1.first)
  {
    vx.push_back(intersect1.second.first);
    vy.push_back(intersect1.second.second);
  }
  pair<bool, pair<double, double>> intersect2 = lines_intersect(ax, ay, bx, by, fx, fy, dx, dy);
  //  bool i2 = lines_intersect(ax, ay, bx, by, fx, fy, dx, dy, &rx, &ry);
  if (intersect2.first)
  {
    vx.push_back(intersect2.second.first);
    vy.push_back(intersect2.second.second);
  }
  pair<bool, pair<double, double>> intersect3 = lines_intersect(ax, ay, bx, by, ex, ey, dx, dy);
  //  bool i3 = lines_intersect(ax, ay, bx, by, ex, ey, dx, dy, &rx, &ry);
  if (intersect3.first)
  {
    vx.push_back(intersect3.second.first);
    vy.push_back(intersect3.second.second);
  }
  pair<bool, pair<double, double>> intersect4 = lines_intersect(ax, ay, bx, by, cx, cy, ex, ey);
  //  bool i4 = lines_intersect(ax, ay, bx, by, cx, cy, ex, ey, &rx, &ry);
  if (intersect4.first)
  {
    vx.push_back(intersect4.second.first);
    vy.push_back(intersect4.second.second);
  }

  //cout << "Rectangle intersections: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
  //cout << "Number of intersections = " << vx.size() << endl;

  double rr = 0.;
  if (vx.size() == 2)
  {
    rr = sqrt((vx[0] - vx[1]) * (vx[0] - vx[1]) + (vy[0] - vy[1]) * (vy[0] - vy[1]));
    //  cout << "Length of intersection = " << *rr << endl;
  }
  if (vx.size() == 1)
  {
    // find which point (A or B) is within the rectangle
    if (ax > cx && ay > cy && ax < dx && ay < dy)  // point A is inside the rectangle
    {
      //cout << "Point A is inside the rectangle." << endl;
      rr = sqrt((vx[0] - ax) * (vx[0] - ax) + (vy[0] - ay) * (vy[0] - ay));
    }
    if (bx > cx && by > cy && bx < dx && by < dy)  // point B is inside the rectangle
    {
      //cout << "Point B is inside the rectangle." << endl;
      rr = sqrt((vx[0] - bx) * (vx[0] - bx) + (vy[0] - by) * (vy[0] - by));
    }
  }

  if (intersect1.first || intersect2.first || intersect3.first || intersect4.first)
  {
    return make_pair(true, rr);
  }
  return make_pair(false, NAN);
}

double PHG4Utils::sA(double r, double x, double y)
{
  // Uses analytic formula for the integral of a circle between limits set by the corner of a rectangle
  // It is called repeatedly to find the overlap area between the circle and rectangle
  // I found this code implementing the integral on a web forum called "ars technica",
  // https://arstechnica.com/civis/viewtopic.php?t=306492
  // posted by "memp"

  double a;

  if (x < 0)
  {
    return -sA(r, -x, y);
  }

  if (y < 0)
  {
    return -sA(r, x, -y);
  }

  if (x > r)
  {
    x = r;
  }

  if (y > r)
  {
    y = r;
  }

  if (x * x + y * y > r * r)
  {
    a = r * r * asin(x / r) + x * sqrt(r * r - x * x) + r * r * asin(y / r) + y * sqrt(r * r - y * y) - r * r * M_PI_2;

    a *= 0.5;
  }
  else
  {
    a = x * y;
  }

  return a;
}

double PHG4Utils::circle_rectangle_intersection(double x1, double y1, double x2, double y2, double mx, double my, double r)
{
  // Find the area of overlap of a circle and rectangle
  // Calls sA, which uses an analytic formula to determine the integral of the circle between limits set by the corners of the rectangle

  // move the rectangle to the frame where the circle is at (0,0)
  x1 -= mx;
  x2 -= mx;
  y1 -= my;
  y2 -= my;

  // {
  //   cout << " mx " << mx << " my " << my << " r " << r << " x1 " << x1 << " x2 " << x2 << " y1 " << y1 << " y2 " << y2 << endl;
  //   cout << " sA21 " << sA(r, x2, y1)
  //        << " sA11 " << sA(r, x1, y1)
  //        << " sA22 " << sA(r, x2, y2)
  //        << " sA12 " << sA(r, x1, y2)
  //        << endl;
  // }

  return sA(r, x2, y1) - sA(r, x1, y1) - sA(r, x2, y2) + sA(r, x1, y2);
}

string PHG4Utils::md5sum(const std::string& filename)
{
  ifstream myfile;
  myfile.open(filename);
  if (!myfile.is_open())
  {
    cout << "Error opening " << filename << endl;
    gSystem->Exit(1);
    exit(1);
  }
  boost::uuids::detail::md5 hash;
  boost::uuids::detail::md5::digest_type digest;
  char c;
  while (myfile.get(c))
  {
    hash.process_bytes(&c, 1);
  }
  myfile.close();
  hash.get_digest(digest);
  const auto charDigest = reinterpret_cast<const char*>(&digest);
  std::string result;
  boost::algorithm::hex(charDigest, charDigest + sizeof(boost::uuids::detail::md5::digest_type), std::back_inserter(result));
  boost::algorithm::to_lower(result);
  return result;
}
