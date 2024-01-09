#ifndef JETBASE_SUBJETINDICESV1_H  
#define JETBASE_SUBJETINDICESV1_H  

#include "SubjetIndices.h"
#include <vector>

class SubjetIndicesv1 : public SubjetIndices
{
  public:

  void Reset() override;

  unsigned int nsubjets() const override { return m_nsubjets; };
  unsigned int index_begin(unsigned int which_jet=0) const override;
  unsigned int index_end  (unsigned int which_jet=0) const override;

  void add_index_pair (unsigned int _begin, unsigned int _end) override;

  private:
    std::vector<unsigned int> v_begin {};
    std::vector<unsigned int> v_end   {};
    unsigned int m_nsubjets {0};
};

#endif // JETBASE_SUBJETINDICESV1_H
