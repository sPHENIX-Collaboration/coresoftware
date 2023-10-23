#include "MbdOut.h"
#include "MbdReturnCodes.h"

#include <cmath>
#include <iostream>

void MbdOut::identify(std::ostream& os) const
{
  os << "virtual MbdOut object" << std::endl;
  return;
}

void MbdOut::Reset()
{
  std::cout << "ERROR MbdOut: Reset() not implemented by daughter class" << std::endl;
  return;
}

int MbdOut::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

Float_t MbdOut::get_zvtx() const
{
  virtual_warning("get_zvtx()");
  return NAN;
}

Float_t MbdOut::get_zvtxerr() const
{
  virtual_warning("get_zvtxerr()");
  return NAN;
}

Float_t MbdOut::get_t0() const
{
  virtual_warning("get_t0()");
  return NAN;
}

//__________________________________________
Float_t MbdOut::get_t0err() const
{
  virtual_warning("get_t0err()");
  return NAN;
}

//__________________________________________
void MbdOut::set_t0(const Float_t /*unused*/, const Float_t /*unused*/)
{
  virtual_warning("set_t0(const Float_t t0, const Float_t t0err)");
  return;
}

//__________________________________________
void MbdOut::set_zvtx(const Float_t /*unused*/, const Float_t /*unused*/)
{
  virtual_warning("set_zvtx(const Float_t vtx, const Float_t vtxerr)");
  return;
}

//__________________________________________
void MbdOut::set_zvtxerr(const Float_t /*unused*/)
{
  virtual_warning("set_zvtxerr(const Float_t vtxerr)");
  return;
}

//________________________________________________________________
void MbdOut::set_arm(const int /*iarm*/, const Short_t /*npmt*/, const Float_t /*energy*/, const Float_t /*timing*/)
{
  virtual_warning("set_arm(const int iMBD, const Short_t npmt, const Float_t energy, const Float_t timing)");
  return;
}

void MbdOut::set_clocks(const Int_t /*evt*/, const UShort_t /*iclk*/, const UShort_t /*ifemclk*/)
{
  virtual_warning("set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk");
  return;
}

Short_t MbdOut::get_npmt(const int /*iarm*/) const
{
  virtual_warning("get_npmt(const int iarm)");
  return MbdReturnCodes::MBD_INVALID_SHORT;
}

Float_t MbdOut::get_q(const int /*iarm*/) const
{
  virtual_warning("get_q(const int iarm)");
  return NAN;
}

Float_t MbdOut::get_time(const int /*iarm*/) const
{
  virtual_warning("get_time(const int iarm)");
  return NAN;
}

Int_t MbdOut::get_evt() const
{
  virtual_warning("get_evt()");
  return 0;
}

UShort_t MbdOut::get_clock() const
{
  virtual_warning("get_clock()");
  return 0;
}

UShort_t MbdOut::get_femclock() const
{
  virtual_warning("get_femclock(const int iarm)");
  return 0;
}

void MbdOut::virtual_warning(const std::string& funcsname) const
{
  std::cout << "MbdOut::" << funcsname << " is virtual, doing nothing" << std::endl;
  return;
}

void MbdOut::FillFromClass(const MbdOut& old)
{
  for (int iarm = 0; iarm < 2; iarm++)
  {
    set_arm(iarm, old.get_npmt(iarm), old.get_q(iarm), old.get_time(iarm));
  }

  set_t0zvtx(old.get_t0(), old.get_t0err(), old.get_zvtx(), old.get_zvtxerr());
}
