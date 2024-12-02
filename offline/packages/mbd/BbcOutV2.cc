#include "BbcOutV2.h"
#include "BbcReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

//______________________________________
BbcOutV2::BbcOutV2() = default;

//______________________________________
void BbcOutV2::Reset()
{
  bz = std::numeric_limits<Float_t>::quiet_NaN();
  bzerr = std::numeric_limits<Float_t>::quiet_NaN();
  bt0 = std::numeric_limits<Float_t>::quiet_NaN();
  bt0err = std::numeric_limits<Float_t>::quiet_NaN();
  bns = 0.;
  bnn = 0.;
  bqs = 0.;
  bqn = 0.;
  bts = std::numeric_limits<Float_t>::quiet_NaN();
  btn = std::numeric_limits<Float_t>::quiet_NaN();
  evt = -1;
  clk = 0;
  femclk = 0;
}

//______________________________________
BbcOutV2::~BbcOutV2() = default;

//______________________________________
int BbcOutV2::isValid() const
{
  // compatible with old invalid setting of -9999.9
  return ((std::isfinite(bt0) && (bt0 > -9999.)) ? 1 : 0);
}

/*
//______________________________________
void BbcOutV2::Reset()
{
  Init();
}
*/

//______________________________________
void BbcOutV2::identify(std::ostream &out) const
{
  out << "identify yourself: I am a BbcOutV2 object" << std::endl;
  out << "Vertex: " << bz << " Error: " << bzerr << std::endl;
  out << "T0: " << bt0 << " Error: " << bt0err << std::endl;
}

//______________________________________
void BbcOutV2::set_t0(const Float_t t0, const Float_t t0err)
{
  bt0 = t0;
  bt0err = t0err;
}

//______________________________________
void BbcOutV2::set_zvtx(const Float_t vtx, const Float_t vtxerr)
{
  bz = vtx;
  bzerr = vtxerr;
}

//______________________________________
void BbcOutV2::set_zvtxerr(const Float_t vtxerr)
{
  bzerr = vtxerr;
}

//______________________________________
void BbcOutV2::set_arm(const int iarm, const Short_t npmt, const Float_t charge, const Float_t timing)
{
  if (iarm == 0)
  {
    bns = npmt;
    bqs = charge;
    bts = timing;
  }
  else if (iarm == 1)
  {
    bnn = npmt;
    bqn = charge;
    btn = timing;
  }
  else
  {
    std::cerr << "BbcOutV2::set_arm(): ERROR, invalid arm " << iarm << std::endl;
  }
}

//______________________________________
void BbcOutV2::set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk)
{
  evt = ievt;
  clk = iclk;
  femclk = ifemclk;
}

//______________________________________
Short_t BbcOutV2::get_npmt(const int iarm) const
{
  return (iarm == 0) ? bns : bnn;
}

//______________________________________
Float_t BbcOutV2::get_q(const int iarm) const
{
  return (iarm == 0) ? bqs : bqn;
}

Float_t BbcOutV2::get_time(const int iarm) const
{
  return (iarm == 0) ? bts : btn;
}

Int_t BbcOutV2::get_evt() const
{
  return evt;
}

UShort_t BbcOutV2::get_clock() const
{
  return clk;
}

UShort_t BbcOutV2::get_femclock() const
{
  return femclk;
}
