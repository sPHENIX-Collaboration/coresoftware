#include "SubjetIndicesv1.h"
#include <phool/phool.h>  // for PHWHERE
#include <iostream>

void SubjetIndicesv1::Reset()
{
  v_begin.clear();
  v_end.clear();
  m_nsubjets = 0;
}

void SubjetIndicesv1::add_index_pair(unsigned int _begin, unsigned int _end)
{
  v_begin.push_back(_begin);
  v_end.push_back(_end);
  m_nsubjets++;
}

unsigned int SubjetIndicesv1::index_begin(unsigned int which_jet) const
{
  if (which_jet >= m_nsubjets)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << " Error -- asking for index of jet indexed(" << which_jet << ") when there are only " << m_nsubjets << std::endl
              << " Returning zero. " << std::endl;
  }
  return v_begin[which_jet];
}

unsigned int SubjetIndicesv1::index_end(unsigned int which_jet) const
{
  if (which_jet >= m_nsubjets)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << " Error -- asking for index of jet indexed(" << which_jet << ") when there are only " << m_nsubjets << std::endl
              << " Returning zero. " << std::endl;
  }
  return v_end[which_jet];
}
