#include "BbcReturnCodes.h"
#include "BbcOutV1.h"
#include "BbcNorthSouthV1.h"

#include <TClonesArray.h>
#include <iostream>

//static const int NPMTBBCV1 = 128;
static const int NBBC = 2;

using namespace std;

ClassImp(BbcOutV1)

//______________________________________
BbcOutV1::BbcOutV1()
{
  Init();
  BbcNS = new TClonesArray("BbcNorthSouthV1",NBBC);
}

//______________________________________
void BbcOutV1::Init()
{
  Bbc_ZVertex   = -99999.9;
  Bbc_dZVertex  = -99999.9;
  Bbc_TimeZero  = -99999.9;
  Bbc_dTimeZero = -99999.9;
}

//______________________________________
BbcOutV1::~BbcOutV1()
{  
  if (BbcNS)
  {
    delete BbcNS;
  }
}

//______________________________________
int BbcOutV1::isValid() const
{
  return((Bbc_TimeZero >-9999.) ? 1 : 0);
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
void BbcOutV1::identify(ostream& out) const
{
  out << "identify yourself: I am a BbcOutV1 object" << endl;
  out << "Vertex: " << Bbc_ZVertex << " Error: " << Bbc_dZVertex << endl;
  out << "T0: " << Bbc_TimeZero << " Error: " << Bbc_dTimeZero << endl;
}

//______________________________________
void BbcOutV1::set_TimeZero(const Float_t t0, const Float_t t0err )
{
  Bbc_TimeZero  = t0;
  Bbc_dTimeZero = t0err;
}

//______________________________________
void BbcOutV1::set_Vertex( const Float_t vtx, const Float_t vtxerr)
{
  Bbc_ZVertex   = vtx;
  Bbc_dZVertex  = vtxerr;
}

//______________________________________
void BbcOutV1::set_dZVertex( const Float_t vtxerr)
{
  Bbc_dZVertex = vtxerr;
}

//______________________________________
void BbcOutV1::AddBbcNS(const int iBBC, const Short_t npmt, const Float_t energy, const Float_t timing)
{
  TClonesArray &bbcns = *BbcNS;
  new(bbcns[iBBC]) BbcNorthSouthV1(npmt, energy,timing);
}

//______________________________________
Short_t BbcOutV1::get_nPMT(const int nBbc) const
{
  BbcNorthSouthV1 *bbcns = (BbcNorthSouthV1*) GetBbcNS()->UncheckedAt(nBbc);
  //  if bbcns=nil (does not exist) return BBC_INVALID_SHORT, else nPMT
  return((bbcns) ? bbcns->get_nPMT() : BBC_INVALID_SHORT);
}

//______________________________________
Float_t BbcOutV1::get_nCharge(const int nBbc) const
{
  BbcNorthSouth *bbcns = (BbcNorthSouthV1*) GetBbcNS()->UncheckedAt(nBbc);
  //  if bbcns=nil (does not exist) return BBC_INVALID_FLOAT, else Energy
  return((bbcns) ? bbcns->get_nCharge() : BBC_INVALID_FLOAT);
}

Float_t BbcOutV1::get_Timing(const int nBbc) const
{
  BbcNorthSouth *bbcns = (BbcNorthSouthV1*) GetBbcNS()->UncheckedAt(nBbc);
  //  if bbcns=nil (does not exist) return BBC_INVALID_FLOAT, else Timing
  return((bbcns) ? bbcns->get_MeanTime() : BBC_INVALID_FLOAT);
}

