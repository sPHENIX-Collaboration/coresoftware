#ifndef PHOOL_PHPOINTERLIST_H
#define PHOOL_PHPOINTERLIST_H

//  Purpose: a template list of pointers
//
//  Description:
//       - The items are held internally as an array of type T.
//       - The list is initialized with a maximum size of 2
//         per default, or with the number given as constructor
//         argument.
//       - If this size is exceeded in the append() function,
//         the array is allocated anew with double size and the
//         old list is copied into this new one.
//       - clear() sets the number of items to zero.
//       - clearAndDestroy() does the same AND deletes all items.
//       - The [] operator always performs a bound-check.
//       - removeLast() returns the pointer to the last item in the
//         array and decrements size by one.
//       - removeAt(i) returns the pointer at position i and rearranges
//         the internal list. This can cost PERFORMANCE in your application.
//       - The output operator '<<' is overloaded for this class.
//         Therefore class T for which the PHPointerList is insantiated must
//         also have an overloaded output operator.
//
//  Author: Matthias Messer

#include <iostream>

template <class T>
class PHPointerList
{
 public:
  explicit PHPointerList(size_t = 2);
  PHPointerList(const PHPointerList<T>&);
  PHPointerList<T>& operator=(const PHPointerList<T>&);
  virtual ~PHPointerList();

 public:
  T* operator[](size_t) const;
  void clear();
  void clearAndDestroy();
  size_t length() const;
  T* removeLast();
  T* removeAt(size_t);
  bool append(T*);
  bool insertAt(T*, size_t);

 private:
  bool grow(size_t = 0);

 private:
  T** items;
  size_t maxNItems;
  size_t nItems;
};

// Implementation of member functions
template <class T>
PHPointerList<T>::PHPointerList(size_t initialSize)
{
  maxNItems = initialSize;
  items = new T*[maxNItems];
  nItems = 0;
}

template <class T>
PHPointerList<T>::PHPointerList(const PHPointerList<T>& l)
{
  *this = l;
}

template <class T>
PHPointerList<T>&
PHPointerList<T>::operator=(const PHPointerList<T>& l)
{
  if (this != &l)
  {
    maxNItems = l.maxNItems;
    grow(l.maxNItems);
    nItems = l.length();
    for (size_t i = 0; i < nItems; ++i)
    {
      items[i] = l[i];
    }
  }
  return *this;
}

template <class T>
PHPointerList<T>::~PHPointerList()
{
  // This deletes the internal list of pointers and NOT the actual objects.
  delete[] items;
}

template <class T>
bool PHPointerList<T>::grow(size_t newSize)
{
  if (newSize == 0)
  {
    newSize = maxNItems * 2;
  }
  T** buffer = items;
  items = new T*[newSize];
  if (items)
  {
    for (size_t i = 0; i < maxNItems; ++i)
    {
      items[i] = buffer[i];
    }
    delete[] buffer;
    maxNItems = newSize;
  }
  else
  {
    std::cout << "PHPointerList<T>::grow: Out of memory?" << std::endl;
    return false;
  }

  return true;
}

template <class T>
inline T* PHPointerList<T>::operator[](size_t i) const
{
  if (i < nItems)
  {
    return items[i];
  }
  else
  {
    std::cout << "PHPointerList<T>::operator[]: nItems exceeded" << std::endl;
    return 0;
  }
}

template <class T>
inline bool
PHPointerList<T>::append(T* item)
{
  if (nItems < maxNItems)
  {
    items[nItems] = item;
    ++nItems;
    return true;
  }
  else
  {
    if (grow())
    {
      items[nItems] = item;
      ++nItems;
      return true;
    }
    else
    {
      std::cout << "PHPointerList<T>::append: max nItems exceeded" << std::endl;
      return false;
    }
  }
}

template <class T>
inline bool
PHPointerList<T>::insertAt(T* item, size_t pos)
{
  // This function inserts item at pos in the internal list
  if (pos > nItems)
  {
    std::cout << "PHPointerList<T>::insertAt: insert beyond nItems" << std::endl;
    return false;
  }

  // Append is used here as a convenient way to let the list grow, if necessary.
  append(item);

  // Now all items are shifted upwards in the list by one, starting at pos.
  for (size_t i = nItems; i > pos; --i)
  {
    items[i] = items[i - 1];
  }

  items[pos] = item;

  return true;
}

template <class T>
inline void
PHPointerList<T>::clear()
{
  nItems = 0;
  items[nItems] = 0;
}

template <class T>
inline void
PHPointerList<T>::clearAndDestroy()
{
  for (size_t i = 0; i < nItems; ++i)
  {
    delete items[i];
  }
  nItems = 0;
  items[nItems] = 0;
}

template <class T>
inline size_t
PHPointerList<T>::length() const
{
  return nItems;
}

template <class T>
inline T*
PHPointerList<T>::removeLast()
{
  if (nItems > 0)
  {
    return items[nItems--];
  }
  else
  {
    std::cout << "PHPointerList<T>::removeLast: no items in list" << std::endl;
    return 0;
  }
}

template <class T>
inline T*
PHPointerList<T>::removeAt(size_t i)
{
  if (i > nItems)
  {
    return 0;
  }

  T* item = items[i];

  for (size_t j = i; j < nItems - 1; ++j)
  {
    items[j] = items[j + 1];
  }
  --nItems;

  return item;
}

// Implementation of external functions.
template <class T>
std::ostream&
operator<<(std::ostream& stream, const PHPointerList<T>& thislist)
{
  for (size_t i = 0; i < thislist.length(); ++i)
  {
    stream << *(thislist[i]) << std::endl;
  }

  return stream;
}

#endif
