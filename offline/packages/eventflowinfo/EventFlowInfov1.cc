#include "EventFlowInfo.h"
#include "EventFlowInfov1.h"
#include <limits>

void EventFlowInfov1::identify(std::ostream& os) const
{
  os << "---------EventFlowInfov1------------------" << std::endl;
  for (const auto& psi : _psi_map) {
    os << "Psi(" << psi.first << ") = " << psi.second << std::endl;
  }
  for (const auto& vn : _vn_map) {
    os << "Vn(" << vn.first << ") = " << vn.second << std::endl;
  }
  return;
}

double EventFlowInfov1::GetPsi(int order) const
{
  auto it = _psi_map.find(order);
  if (it == _psi_map.end()) 
  {
    std::cerr << "EventFlowInfov1::GetPsi: order " << order << " not found in psi_map" << std::endl;
    return NAN;
  }
  return it->second;
}


double EventFlowInfov1::GetVn(int order) const
{
  auto it = _vn_map.find(order);
  if (it == _vn_map.end()) 
  {
    std::cerr << "EventFlowInfov1::GetVn: order " << order << " not found in vn_map" << std::endl;
    return NAN;
  }
  return it->second;

}