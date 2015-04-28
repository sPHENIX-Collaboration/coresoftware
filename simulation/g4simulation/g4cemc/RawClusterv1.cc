#include "RawClusterv1.h"

using namespace std;

ClassImp(RawClusterv1)

RawClusterv1::RawClusterv1() : RawCluster(), _eta(0.0), _phi(0.0), _energy(0.0)
{}

void
RawClusterv1::Reset()
{
  _eta = 0.0;
  _phi = 0.0;
  _energy = 0.0;
  _towers.clear();
}

std::pair<int,int>
RawClusterv1::getTowerBin(const unsigned int itower) const
{
  if (itower < _towers.size())
    {
      return  _towers[itower];
    }
  return make_pair(-9999,-9999);
}

