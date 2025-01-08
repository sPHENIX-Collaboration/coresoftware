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

Float_t BbcOut::get_zvtx() const
{
  virtual_warning("get_zvtx()");
  return NAN;
}

Float_t BbcOut::get_zvtxerr() const
{
  virtual_warning("get_zvtxerr()");
  return NAN;
}

Float_t BbcOut::get_t0() const
{
  virtual_warning("get_t0()");
  return NAN;
}

//__________________________________________
Float_t BbcOut::get_t0err() const
{
  virtual_warning("get_t0err()");
  return NAN;
}

//__________________________________________
void BbcOut::set_t0(const Float_t /*unused*/, const Float_t /*unused*/)
{
  virtual_warning("set_t0(const Float_t t0, const Float_t t0err)");
  return;
}

//__________________________________________
void BbcOut::set_zvtx(const Float_t /*unused*/, const Float_t /*unused*/)
{
  virtual_warning("set_zvtx(const Float_t vtx, const Float_t vtxerr)");
  return;
}

//__________________________________________
void BbcOut::set_zvtxerr(const Float_t /*unused*/)
{
  virtual_warning("set_zvtxerr(const Float_t vtxerr)");
  return;
}

//________________________________________________________________
void BbcOut::set_arm(const int /*iarm*/, const Short_t /*npmt*/, const Float_t /*energy*/, const Float_t /*timing*/)
{
  virtual_warning("set_arm(const int iBBC, const Short_t npmt, const Float_t energy, const Float_t timing)");
  return;
}

void BbcOut::set_clocks(const Int_t /*evt*/, const UShort_t /*iclk*/, const UShort_t /*ifemclk*/)
{
  virtual_warning("set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk");
  return;
}

Short_t BbcOut::get_npmt(const int /*iarm*/) const
{
  virtual_warning("get_npmt(const int iarm)");
  return BbcReturnCodes::BBC_INVALID_SHORT;
}

Float_t BbcOut::get_q(const int /*iarm*/) const
{
  virtual_warning("get_q(const int iarm)");
  return NAN;
}

Float_t BbcOut::get_time(const int /*iarm*/) const
{
  virtual_warning("get_time(const int iarm)");
  return NAN;
}

Int_t BbcOut::get_evt() const
{
  virtual_warning("get_evt()");
  return 0;
}

UShort_t BbcOut::get_clock() const
{
  virtual_warning("get_clock()");
  return 0;
}

UShort_t BbcOut::get_femclock() const
{
  virtual_warning("get_femclock(const int iarm)");
  return 0;
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
