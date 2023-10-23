#include "MbdOutV1.h"
#include "MbdReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

//______________________________________
MbdOutV1::MbdOutV1()
{
}

//______________________________________
void MbdOutV1::Reset()
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
}

//______________________________________
MbdOutV1::~MbdOutV1()
{
}

//______________________________________
int MbdOutV1::isValid() const
{
// compatible with old invalid setting of -9999.9
  return ((std::isfinite(bt0) && (bt0 > -9999.)) ? 1 : 0);
}

//______________________________________
void MbdOutV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdOutV1 object" << std::endl;
  out << "Vertex: " << bz << " Error: " << bzerr << std::endl;
  out << "T0: " << bt0 << " Error: " << bt0err << std::endl;
}

//______________________________________
void MbdOutV1::set_t0(const Float_t t0, const Float_t t0err)
{
  bt0 = t0;
  bt0err = t0err;
}

//______________________________________
void MbdOutV1::set_zvtx(const Float_t vtx, const Float_t vtxerr)
{
  bz = vtx;
  bzerr = vtxerr;
}

//______________________________________
void MbdOutV1::set_zvtxerr(const Float_t vtxerr)
{
  bzerr = vtxerr;
}

//______________________________________
void MbdOutV1::set_arm(const int iarm, const Short_t npmt, const Float_t charge, const Float_t timing)
{
  if (iarm==0)
  {
    bns = npmt;
    bqs = charge;
    bts = timing;
  }
  else if (iarm==1)
  {
    bnn = npmt;
    bqn = charge;
    btn = timing;
  }
  else
  {
    std::cerr << "MbdOutV1::set_arm(): ERROR, invalid arm " << iarm << std::endl;
  }
}

//______________________________________
Short_t MbdOutV1::get_npmt(const int iarm) const
{
  return (iarm==0) ? bns : bnn;
}

//______________________________________
Float_t MbdOutV1::get_q(const int iarm) const
{
  return (iarm==0) ? bqs : bqn;
}

Float_t MbdOutV1::get_time(const int iarm) const
{
  return (iarm==0) ? bts : btn;
}

