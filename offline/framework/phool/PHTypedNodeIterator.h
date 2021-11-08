#ifndef PHOOL_PHTYPEDNODEITERATOR_H
#define PHOOL_PHTYPEDNODEITERATOR_H

#include "PHNodeIterator.h"

#include <cstddef>

class TObject;
template <typename T>
class PHIODataNode;
class PHCompositeNode;

//  A special PHOOL node iterator that simplifies finding and adding
//  data nodes.  This version is a templated class, to allow its use
//  interactively in CINT.  (CINT only allows function templates in
//  very recent versions, which we aren't using yet.)
//  @author Kyle Pope

template <class T>
class PHTypedNodeIterator : public PHNodeIterator
{
 public:
  /// Constructor
  explicit PHTypedNodeIterator(PHCompositeNode* n)
    : PHNodeIterator(n)
  {
    myIODataNode = nullptr;
  }

  PHTypedNodeIterator()
    : PHNodeIterator()
  {
    myIODataNode = nullptr;
  }

  T& operator*();

  /// Destructor
  //  virtual ~PHTypedNodeIterator();  // Need a virtual ~PHNodeIterator() !
  ~PHTypedNodeIterator() override {}
  /**
   * Finds an IODataNode of name "name" containing data of type "T".
   * A null pointer will be returned if the node is not found, or if
   * it contains data of the wrong type.
   */
  PHIODataNode<T>* FindIODataNode(const char* name);

  PHIODataNode<T>* find(const char* name);

  /**
   * Adds a data node called "name" to the tree, and inserts "data".
   * The data node is added at the current "directory" of this iterator
   * object, so remember to "cd" to the desired location in the tree!
   *
   */
  bool AddIODataNode(T* data, const char* name);

  bool insert(T* data, const char* name);

 protected:
  PHIODataNode<T>* myIODataNode;
};

template <class T>
T&
    PHTypedNodeIterator<T>::operator*()
{
  return *(myIODataNode->getData());
}

template <class T>
PHIODataNode<T>*
PHTypedNodeIterator<T>::FindIODataNode(const char* name)
{
  return find(name);
}

template <class T>
PHIODataNode<T>*
PHTypedNodeIterator<T>::find(const char* name)
{
  // TODO:  also check that "name" is not a null string!
  if (!name)
  {
    return 0;
  }

  // Can't do dynamic_cast here; it fails if node was created as
  // PHIODataNode<X> instead of PHIODataNode<T>, even if T is a
  // derived class of X!  In general, T -> X does not imply A<T> ->
  // A<X>.  ("->" denotes "derives from", and "A" is any template
  // class)
  myIODataNode =
      static_cast<PHIODataNode<T>*>(findFirst("PHIODataNode", name));

  return myIODataNode;
}

template <class T>
bool PHTypedNodeIterator<T>::AddIODataNode(T* data, const char* name)
{
  return insert(data, name);
}

template <class T>
bool PHTypedNodeIterator<T>::insert(T* data, const char* name)
{
  // TODO:  also check that "name" is not a null string!
  if (!name)
  {
    return false;
  }

  // For IODataNode, ought to check (if possible) that T derives from
  // TObject.  Will typeid() give us this info?
  PHIODataNode<T>* n = new PHIODataNode<T>(data, name);
  if (!n)
  {
    // problem creating node?
    return false;
  }

  return addNode(n);
}

// Typedef to simplify notation.
typedef PHTypedNodeIterator<TObject> PHRootNodeIterator;

#endif
