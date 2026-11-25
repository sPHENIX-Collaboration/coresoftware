#include "TpcClusterdxMapv1.h"

#include <map>

TpcClusterdxMapv1::TpcClusterdxMapv1()
{
}

TpcClusterdxMapv1::~TpcClusterdxMapv1()
{
  Cluster_dx.clear();
}

void TpcClusterdxMapv1::Reset()
{
  Cluster_dx.clear();
}

void TpcClusterdxMapv1::identify(std::ostream &os) const
{
  os << "TpcClusterdxMapv1" << std::endl;
  os << "containing " << Cluster_dx.size() << " dx values." << std::endl;
}

int TpcClusterdxMapv1::isValid() const
{
  return Cluster_dx.size();
}

double TpcClusterdxMapv1::Add_Cluster_dx(uint64_t index, double cluster_dx)
{
  Cluster_dx.insert({index,cluster_dx});
  return cluster_dx;
}

double TpcClusterdxMapv1::get_Cluster_dx(uint64_t index)
{
  return (double) Cluster_dx.at(index);
}

std::map<uint64_t,double> TpcClusterdxMapv1::get_map()
{
  return Cluster_dx;
}
