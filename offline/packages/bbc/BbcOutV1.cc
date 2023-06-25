#include "BbcOutV1.h"
#include "BbcNorthSouthV1.h"
#include "BbcReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

static const int NBBC = 2;

//______________________________________
BbcOutV1::BbcOutV1()
{
  Init();
  BbcNS = new TClonesArray("BbcNorthSouthV1", NBBC);
}

//______________________________________
void BbcOutV1::Init()
{
  Bbc_ZVertex = std::numeric_limits<float>::quiet_NaN();
  Bbc_dZVertex = std::numeric_limits<float>::quiet_NaN();
  Bbc_TimeZero = std::numeric_limits<float>::quiet_NaN();
  Bbc_dTimeZero = std::numeric_limits<float>::quiet_NaN();
}

//______________________________________
BbcOutV1::~BbcOutV1()
{
  delete BbcNS;
}

//______________________________________
int BbcOutV1::isValid() const
{
// compatible with old invalid setting of -9999.9
  return ((std::isfinite(Bbc_TimeZero) && (Bbc_TimeZero > -9999.)) ? 1 : 0);
}

//______________________________________
void BbcOutV1::Reset()
{
  Init();
  if (BbcNS)
  {
    BbcNS->Clear();
  }
}

//______________________________________
void BbcOutV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a BbcOutV1 object" << std::endl;
  out << "Vertex: " << Bbc_ZVertex << " Error: " << Bbc_dZVertex << std::endl;
  out << "T0: " << Bbc_TimeZero << " Error: " << Bbc_dTimeZero << std::endl;
}

//______________________________________
void BbcOutV1::set_TimeZero(const float t0, const float t0err)
{
  Bbc_TimeZero = t0;
  Bbc_dTimeZero = t0err;
}

//______________________________________
void BbcOutV1::set_Vertex(const float vtx, const float vtxerr)
{
  Bbc_ZVertex = vtx;
  Bbc_dZVertex = vtxerr;
}

//______________________________________
void BbcOutV1::set_dZVertex(const float vtxerr)
{
  Bbc_dZVertex = vtxerr;
}

//______________________________________
void BbcOutV1::AddBbcNS(const int iBBC, const short npmt, const float energy, const float timing)
{
  TClonesArray &bbcns = *BbcNS;
  new (bbcns[iBBC]) BbcNorthSouthV1(npmt, energy, timing);
}

//______________________________________
short BbcOutV1::get_nPMT(const int nBbc) const
{
  BbcNorthSouthV1 *bbcns = static_cast<BbcNorthSouthV1 *> (GetBbcNS()->UncheckedAt(nBbc));
  //  if bbcns=null (does not exist) return BbcReturnCodes::BBC_INVALID_SHORT, else nPMT
  return ((bbcns) ? bbcns->get_nPMT() : BbcReturnCodes::BBC_INVALID_SHORT);
}

//______________________________________
float BbcOutV1::get_nCharge(const int nBbc) const
{
  BbcNorthSouth *bbcns = static_cast<BbcNorthSouthV1 *> (GetBbcNS()->UncheckedAt(nBbc));
  //  if bbcns=null (does not exist) return BbcReturnCodes::BBC_INVALID_FLOAT, else Energy
  return ((bbcns) ? bbcns->get_nCharge() : BbcReturnCodes::BBC_INVALID_FLOAT);
}

float BbcOutV1::get_Timing(const int nBbc) const
{
  BbcNorthSouth *bbcns = static_cast<BbcNorthSouthV1 *> (GetBbcNS()->UncheckedAt(nBbc));
  //  if bbcns=null (does not exist) return BbcReturnCodes::BBC_INVALID_FLOAT, else Timing
  return ((bbcns) ? bbcns->get_MeanTime() : BbcReturnCodes::BBC_INVALID_FLOAT);
}
