#ifndef __PHPOINTERLISTITERATOR_H__
#define __PHPOINTERLISTITERATOR_H__

//  Declaration of class PHPointerListIterator
//  Purpose: iterator for access to a PHPointerList
//  Author: Matthias Messer

#include "PHPointerList.h"

template <class T>
class PHPointerListIterator
{
 public:
  PHPointerListIterator(const PHPointerList<T>&);
  virtual ~PHPointerListIterator() {}
 public:
  T* operator()();
  void operator--();
  void reset();
  size_t pos() const { return index; }
 protected:
  PHPointerListIterator()
    : list(0)
    , index(0)
  {
  }

 private:
  const PHPointerList<T>& list;
  size_t index;
};

template <class T>
PHPointerListIterator<T>::PHPointerListIterator(const PHPointerList<T>& lis)
  : list(lis)
{
  reset();
}

template <class T>
T* PHPointerListIterator<T>::operator()()
{
  index++;
  if (index < list.length())
  {
    return list[index];
  }
  else
  {
    return 0;
  }
}

template <class T>
void
    PHPointerListIterator<T>::operator--()
{
  --index;
}

template <class T>
void PHPointerListIterator<T>::reset()
{
  index = ~(size_t) 0;
}

#endif /* __PHPOINTERLISTITERATOR_H__ */
