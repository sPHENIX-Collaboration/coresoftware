#include "BbcOut.h"
#include "BbcReturnCodes.h"

#include <cmath>
#include <iostream>

void BbcOut::identify(std::ostream& os) const
{
  os << "virtual BbcOut object" << std::endl;
  return;
}

void BbcOut::Reset()
{
  std::cout << "ERROR BbcOut: Reset() not implemented by daughter class" << std::endl;
  return;
}

int BbcOut::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

float BbcOut::get_VertexPoint() const
{
  virtual_warning("get_VertexPoint()");
  return NAN;
}

float BbcOut::get_dVertexPoint() const
{
  virtual_warning("get_dVertexPoint()");
  return NAN;
}

float BbcOut::get_TimeZero() const
{
  virtual_warning("get_TimeZero()");
  return NAN;
}

//__________________________________________
float BbcOut::get_dTimeZero() const
{
  virtual_warning("get_dTimeZero()");
  return NAN;
}

//__________________________________________
void BbcOut::set_TimeZero(const float /*unused*/, const float /*unused*/)
{
  virtual_warning("set_TimeZero(const float t0, const float t0err)");
  return;
}

//__________________________________________
void BbcOut::set_Vertex(const float /*unused*/, const float /*unused*/)
{
  virtual_warning("set_Vertex(const float vtx, const float vtxerr)");
  return;
}

//__________________________________________
void BbcOut::set_dZVertex(const float /*unused*/)
{
  virtual_warning("set_dZVertex(const float vtxerr)");
  return;
}

//________________________________________________________________
void BbcOut::AddBbcNS(const int /*nBbc*/, const short /*npmt*/, const float /*energy*/, const float /*timing*/)
{
  virtual_warning("AddBbcNS(const int iBBC, const short npmt, const float energy, const float timing)");
  return;
}

short BbcOut::get_nPMT(const int /*nBbc*/) const
{
  virtual_warning("get_nPMT(const int nBbc)");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

float BbcOut::get_nCharge(const int /*nBbc*/) const
{
  virtual_warning("get_nCharge(const int nBbc)");
  return NAN;
}

float BbcOut::get_Timing(const int /*nBbc*/) const
{
  virtual_warning("get_Timing(const int nBbc)");
  return NAN;
}

void BbcOut::virtual_warning(const std::string& funcsname) const
{
  std::cout << "BbcOut::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}

void BbcOut::FillFromClass(const BbcOut& old)
{
  for (int iarm = 0; iarm < 2; iarm++)
  {
    AddBbcNS(iarm, old.get_nPMT(iarm), old.get_nCharge(iarm), old.get_Timing(iarm));
  }

  set_TimeVertex(old.get_TimeZero(), old.get_dTimeZero(), old.get_VertexPoint(), old.get_dVertexPoint());
}
