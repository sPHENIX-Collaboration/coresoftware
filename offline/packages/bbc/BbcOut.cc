#include "BbcReturnCodes.h"
#include "BbcOut.h"

#include <cmath>
#include <iostream>

ClassImp(BbcOut)

void BbcOut::identify(std::ostream& os) const
{
  os << "virtual BbcOut object";
  return ;
}

void BbcOut::Reset()
{
  std::cout << "ERROR BbcOut: Reset() not implemented by daughter class" << std::endl;
  return ;
}

int BbcOut::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

Float_t BbcOut::get_VertexPoint() const
{
  virtual_warning("get_VertexPoint()");
  return NAN;
}

Float_t BbcOut::get_dVertexPoint() const
{
  virtual_warning("get_dVertexPoint()");
  return NAN;
}

Float_t BbcOut::get_TimeZero() const
{
  virtual_warning("get_TimeZero()");
  return NAN;
}

//__________________________________________
Float_t BbcOut::get_dTimeZero() const
{
  virtual_warning("get_dTimeZero()");
  return NAN;
}

//__________________________________________
void BbcOut::set_TimeZero(const Float_t, const Float_t)
{
  virtual_warning("set_TimeZero(const Float_t t0, const Float_t t0err)");
  return ;
}

//__________________________________________
void BbcOut::set_Vertex( const Float_t, const Float_t )
{
  virtual_warning("set_Vertex(const Float_t vtx, const Float_t vtxerr)");
  return ;
}

//__________________________________________
void BbcOut::set_dZVertex(const Float_t )
{
  virtual_warning("set_dZVertex(const Float_t vtxerr)");
  return ;
}

//________________________________________________________________
void BbcOut::AddBbcNS(const int /*nBbc*/, const Short_t /*npmt*/, const Float_t /*energy*/, const Float_t /*timing*/)
{
  virtual_warning("AddBbcNS(const int iBBC, const Short_t npmt, const Float_t energy, const Float_t timing)");
  return ;
}

Short_t BbcOut::get_nPMT(const int /*nBbc*/) const
{
  virtual_warning("get_nPMT(const int nBbc)");
  return BBC_INVALID_SHORT;
}

Float_t BbcOut::get_nCharge(const int /*nBbc*/) const
{
  virtual_warning("get_nCharge(const int nBbc)");
  return NAN;
}

Float_t BbcOut::get_Timing(const int /*nBbc*/) const
{
  virtual_warning("get_Timing(const int nBbc)");
  return NAN;
}

void BbcOut::virtual_warning(const char *funcsname) const
{
  std::cout << "BbcOut::" << funcsname << " is virtual, doing nothing" << std::endl;
  return ;
}

void BbcOut::FillFromClass(const BbcOut& old)
{
  for(int iarm = 0; iarm < 2; iarm++)
  {
    AddBbcNS( iarm, old.get_nPMT(iarm), old.get_nCharge(iarm), old.get_Timing(iarm) );
  }

  set_TimeVertex( old.get_TimeZero(), old.get_dTimeZero(), old.get_VertexPoint(), old.get_dVertexPoint() );

}
