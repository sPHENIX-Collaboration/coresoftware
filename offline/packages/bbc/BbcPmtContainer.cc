#include <iostream>
#include "phool/phool.h"
#include "BbcPmtContainer.h"
#include "BbcReturnCodes.h"

using namespace std;

ClassImp(BbcPmtContainer)

void BbcPmtContainer::identify(ostream& os) const
{
  os << "virtual BbcPmtContainer object";
  return ;
}

void BbcPmtContainer::Reset()
{
  cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << endl;
  return ;
}

int BbcPmtContainer::isValid() const
{
  virtual_warning("isValid()");
  return 0;
}

void BbcPmtContainer::set_npmt(const Short_t /*ival*/)
{
  virtual_warning("set_npmt(const short ival)");
  return ;
}

short BbcPmtContainer::get_npmt() const
{
  virtual_warning("get_npmt()");
  return BBC_INVALID_SHORT;
}

Float_t BbcPmtContainer::get_adc(const int /*iPmt*/) const
{
  virtual_warning("get_Adc(const short iPmt)");
  return BBC_INVALID_SHORT;
}

Float_t BbcPmtContainer::get_tdc0(const int /*iPmt*/) const
{
  virtual_warning("get_Tdc0(const short iPmt)");
  return BBC_INVALID_SHORT;
}

Float_t BbcPmtContainer::get_tdc1(const int /*iPmt*/) const
{
  virtual_warning("get_Tdc1(const short iPmt)");
  return BBC_INVALID_SHORT;
}

void BbcPmtContainer::AddBbcPmt(const Short_t /*ipmt*/, const Float_t /*adc*/, const Float_t /*tdc0*/, const Float_t /*tdc1*/)
{
  virtual_warning("AddBbcPmtHit(const Short_t ipmt, const Short_t adc, const Float_t tdc0, const Float_t tdc1)");
  return ;
}

void BbcPmtContainer::virtual_warning(const char *funcsname) const
{
  cout << "BbcPmtContainer::" << funcsname << " is virtual, doing nothing" << endl;
  return ;
}
