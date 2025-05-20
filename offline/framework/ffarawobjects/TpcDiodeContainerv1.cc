#include "TpcDiodeContainerv1.h"
#include "TpcDiodev1.h"

#include <TClonesArray.h>

static const int NTPCDIODES = 10000;

TpcDiodeContainerv1::TpcDiodeContainerv1()
  : TpcDiodesTCArray(new TClonesArray("TpcDiodev1", NTPCDIODES))
{
}

TpcDiodeContainerv1::~TpcDiodeContainerv1()
{
  TpcDiodesTCArray->Clear("C");
  delete TpcDiodesTCArray;
}

void TpcDiodeContainerv1::Reset()
{
  TpcDiodesTCArray->Clear("C");
  TpcDiodesTCArray->Expand(NTPCDIODES);
}

void TpcDiodeContainerv1::identify(std::ostream &os) const
{
  os << "TpcDiodeContainerv1" << std::endl;
  os << "containing " << TpcDiodesTCArray->GetEntriesFast() << " Tpc diodes" << std::endl;
  // TpcDiode *tpcdiode = static_cast<TpcDiode *>(TpcDiodesTCArray->At(0));
  //  if (tpcdiode)
  //  {
  //    os << "for beam clock: " << std::hex << tpcdiode->get_bco() << std::dec << std::endl;
  //  }
}

int TpcDiodeContainerv1::isValid() const
{
  return TpcDiodesTCArray->GetSize();
}

unsigned int TpcDiodeContainerv1::get_ndiodes()
{
  return TpcDiodesTCArray->GetEntriesFast();
}

TpcDiode *TpcDiodeContainerv1::AddDiode()
{
  TpcDiode *newdiode = new ((*TpcDiodesTCArray)[TpcDiodesTCArray->GetLast() + 1]) TpcDiodev1();
  return newdiode;
}

TpcDiode *TpcDiodeContainerv1::AddDiode(TpcDiode *tpcdiode)
{
  TpcDiode *newdiode = new ((*TpcDiodesTCArray)[TpcDiodesTCArray->GetLast() + 1]) TpcDiodev1(tpcdiode);
  return newdiode;
}

TpcDiode *TpcDiodeContainerv1::get_diode(unsigned int index)
{
  return (TpcDiode *) TpcDiodesTCArray->At(index);
}

unsigned int TpcDiodeContainerv1::get_Laser()
{
  int laser = -1;
  std::vector<int> nlasers;
  int threshold = 200;
  for (int c = 0; c < 32; c++)
  {
    if ((c > 3 && c < 16) || c > 19)
    {
      continue;
    }
    TpcDiode *EMon = get_diode(c);
    int maxadc = EMon->get_maxadc();
    if (maxadc > threshold)
    {
      laser = c;
      nlasers.push_back(c);
    }
  }
  if (nlasers.size() > 1)
  {
    std::cout << "More than one laser fired!" << std::endl;
    for (int nlaser : nlasers)
    {
      std::cout << "Laser " << nlaser << " fired" << std::endl;
    }
    return -1;
  }
  if (laser < 0)
  {
    std::cout << "No laser fired!" << std::endl;
    return laser;
  }
  return laser;
}

// code below this point needs to be rewritten!!!

std::vector<TpcDiode *> TpcDiodeContainerv1::get_PO1()
{
  std::vector<TpcDiode *> PO1;

  unsigned int laser = get_Laser();

  if (laser == 0 || laser == 1 || laser == 2 || laser == 3)
  {
    PO1.push_back(get_diode(4));
    PO1.push_back(get_diode(5));
    PO1.push_back(get_diode(6));
    PO1.push_back(get_diode(7));
  }

  else if (laser == 16 || laser == 17 || laser == 18 || laser == 19)
  {
    PO1.push_back(get_diode(20));
    PO1.push_back(get_diode(21));
    PO1.push_back(get_diode(22));
    PO1.push_back(get_diode(23));
  }
  else
  {
    std::cout << "No laser fired in this event!" << std::endl;
  }
  return PO1;
}

std::vector<TpcDiode *> TpcDiodeContainerv1::get_PO2()
{
  std::vector<TpcDiode *> PO2;

  unsigned int laser = get_Laser();

  if (laser == 0 || laser == 1 || laser == 2 || laser == 3)
  {
    PO2.push_back(get_diode(8));
    PO2.push_back(get_diode(9));
    PO2.push_back(get_diode(10));
    PO2.push_back(get_diode(11));
  }

  else if (laser == 16 || laser == 17 || laser == 18 || laser == 19)
  {
    PO2.push_back(get_diode(24));
    PO2.push_back(get_diode(25));
    PO2.push_back(get_diode(26));
    PO2.push_back(get_diode(27));
  }
  else
  {
    std::cout << "No laser fired in this event!" << std::endl;
  }
  return PO2;
}

std::vector<TpcDiode *> TpcDiodeContainerv1::get_EGG()
{
  std::vector<TpcDiode *> EGG;

  unsigned int laser = get_Laser();

  if (laser == 0 || laser == 1 || laser == 2 || laser == 3)
  {
    EGG.push_back(get_diode(12));
    EGG.push_back(get_diode(13));
    EGG.push_back(get_diode(14));
    EGG.push_back(get_diode(15));
  }

  else if (laser == 16 || laser == 17 || laser == 18 || laser == 19)
  {
    EGG.push_back(get_diode(28));
    EGG.push_back(get_diode(29));
    EGG.push_back(get_diode(30));
    EGG.push_back(get_diode(31));
  }
  else
  {
    std::cout << "No laser fired in this event!" << std::endl;
  }
  return EGG;
}
