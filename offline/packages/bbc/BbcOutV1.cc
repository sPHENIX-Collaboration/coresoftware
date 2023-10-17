#include "BbcOutV1.h"
#include "BbcReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

//______________________________________
BbcOutV1::BbcOutV1()
{
  Init();
}

//______________________________________
void BbcOutV1::Init()
{
  bz = std::numeric_limits<float>::quiet_NaN();
  bzerr = std::numeric_limits<float>::quiet_NaN();
  bt0 = std::numeric_limits<float>::quiet_NaN();
  bt0err = std::numeric_limits<float>::quiet_NaN();
  bns = 0.;
  bnn = 0.;
  bqs = 0.;
  bqn = 0.;
  bts = std::numeric_limits<float>::quiet_NaN();
  btn = std::numeric_limits<float>::quiet_NaN();
}

//______________________________________
BbcOutV1::~BbcOutV1()
{
}

//______________________________________
int BbcOutV1::isValid() const
{
// compatible with old invalid setting of -9999.9
  return ((std::isfinite(bt0) && (bt0 > -9999.)) ? 1 : 0);
}

//______________________________________
void BbcOutV1::Reset()
{
  Init();
}

//______________________________________
void BbcOutV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a BbcOutV1 object" << std::endl;
  out << "Vertex: " << bz << " Error: " << bzerr << std::endl;
  out << "T0: " << bt0 << " Error: " << bt0err << std::endl;
}

//______________________________________
void BbcOutV1::set_t0(const float t0, const float t0err)
{
  bt0 = t0;
  bt0err = t0err;
}

//______________________________________
void BbcOutV1::set_zvtx(const float vtx, const float vtxerr)
{
  bz = vtx;
  bzerr = vtxerr;
}

//______________________________________
void BbcOutV1::set_zvtxerr(const float vtxerr)
{
  bzerr = vtxerr;
}

//______________________________________
void BbcOutV1::set_arm(const int iarm, const short npmt, const float charge, const float timing)
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
    std::cerr << "BbcOutV1::set_arm(): ERROR, invalid arm " << iarm << std::endl;
  }
}

//______________________________________
short BbcOutV1::get_npmt(const int iarm) const
{
  return (iarm==0) ? bns : bnn;
}

//______________________________________
float BbcOutV1::get_q(const int iarm) const
{
  return (iarm==0) ? bqs : bqn;
}

float BbcOutV1::get_time(const int iarm) const
{
  return (iarm==0) ? bts : btn;
}
