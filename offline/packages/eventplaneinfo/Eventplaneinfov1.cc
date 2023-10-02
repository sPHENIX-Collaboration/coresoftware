#include "Eventplaneinfov1.h"

#include <cmath>

Eventplaneinfov1::Eventplaneinfov1(const Eventplaneinfo::EPTYPE id)
  : _id(id)
{
}

void Eventplaneinfov1::identify(std::ostream& os) const
{
  os << "---Eventplaneinfov1-----------------------" << std::endl;

  os << "epid: " << get_id();
 
  os <<"\t raw psi_2 "<< get_psi_raw(2) <<std::endl;
    
  os << " list of ep ids: " << std::endl;
  for (ConstEpIter iter = begin_ep_ids(); iter != end_ep_ids(); ++iter)
  {
    os << "  " << iter->first << " => " << iter->second << std::endl;
  }
  os << "-----------------------------------------------" << std::endl;

  return;
}



