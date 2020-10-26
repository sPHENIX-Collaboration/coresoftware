#ifndef PHOOL_PHPOINTERLISTITERATOR_H
#define PHOOL_PHPOINTERLISTITERATOR_H

//  Declaration of class PHPointerListIterator
//  Purpose: iterator for access to a PHPointerList
//  Author: Matthias Messer

#include "PHPointerList.h"

template <class T>
class PHPointerListIterator
{
 public:
  explicit PHPointerListIterator(const PHPointerList<T>&);
  virtual ~PHPointerListIterator() {}
  T* operator()();
  void operator--();
  void reset();
  size_t pos() const { return m_Index; }

 private:
  PHPointerListIterator() = delete;
  const PHPointerList<T>& m_List;
  size_t m_Index;
};

template <class T>
PHPointerListIterator<T>::PHPointerListIterator(const PHPointerList<T>& lis)
  : m_List(lis)
{
  reset();
}

template <class T>
T* PHPointerListIterator<T>::operator()()
{
  m_Index++;
  if (m_Index < m_List.length())
  {
    return m_List[m_Index];
  }
  else
  {
    return 0;
  }
}

template <class T>
void PHPointerListIterator<T>::operator--()
{
  --m_Index;
}

template <class T>
void PHPointerListIterator<T>::reset()
{
  m_Index = ~(size_t) 0;
}

#endif  // PHOOL_PHPOINTERLISTITERATOR_H
