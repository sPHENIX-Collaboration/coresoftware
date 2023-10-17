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

float BbcOut::get_zvtx() const
{
  virtual_warning("get_zvtx()");
  return NAN;
}

float BbcOut::get_zvtxerr() const
{
  virtual_warning("get_zvtxerr()");
  return NAN;
}

float BbcOut::get_t0() const
{
  virtual_warning("get_t0()");
  return NAN;
}

//__________________________________________
float BbcOut::get_t0err() const
{
  virtual_warning("get_t0err()");
  return NAN;
}

//__________________________________________
void BbcOut::set_t0(const float /*unused*/, const float /*unused*/)
{
  virtual_warning("set_t0(const float t0, const float t0err)");
  return;
}

//__________________________________________
void BbcOut::set_zvtx(const float /*unused*/, const float /*unused*/)
{
  virtual_warning("set_zvtx(const float vtx, const float vtxerr)");
  return;
}

//__________________________________________
void BbcOut::set_zvtxerr(const float /*unused*/)
{
  virtual_warning("set_zvtxerr(const float vtxerr)");
  return;
}

//________________________________________________________________
void BbcOut::set_arm(const int /*iarm*/, const short /*npmt*/, const float /*energy*/, const float /*timing*/)
{
  virtual_warning("set_arm(const int iBBC, const short npmt, const float energy, const float timing)");
  return;
}

short BbcOut::get_npmt(const int /*iarm*/) const
{
  virtual_warning("get_npmt(const int iarm)");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

float BbcOut::get_q(const int /*iarm*/) const
{
  virtual_warning("get_q(const int iarm)");
  return NAN;
}

float BbcOut::get_time(const int /*iarm*/) const
{
  virtual_warning("get_time(const int iarm)");
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
    set_arm(iarm, old.get_npmt(iarm), old.get_q(iarm), old.get_time(iarm));
  }

  set_t0zvtx(old.get_t0(), old.get_t0err(), old.get_zvtx(), old.get_zvtxerr());
}
