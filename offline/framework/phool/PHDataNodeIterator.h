#ifndef PHOOL_PHDATANODEITERATOR_H
#define PHOOL_PHDATANODEITERATOR_H

#include "PHIODataNode.h"
#include "PHNodeIterator.h"

#include <cstddef>

/**
 * A special PHOOL node iterator that simplifies finding and adding data
 * nodes.
 *
 * The methods are templated, rather than the class, to allow the same
 * iterator object to be used with different data node types.
 *
 * @author Kyle Pope
 */

class PHDataNodeIterator : public PHNodeIterator
{
 public:
  /// Constructor
  PHDataNodeIterator(PHCompositeNode* node);

  /// Destructor
  ~PHDataNodeIterator() override {}
  /**
   * Finds an IODataNode of name "name" containing data of type "T".
   * A null pointer will be returned if the node is not found, or if
   * it contains data of the wrong type.  Note that the return
   * variable is also the first argument; this is necessary to resolve
   * the template type.
   */
  template <class T>
  PHIODataNode<T>* FindIODataNode(PHIODataNode<T>* node,
                                  const char* name);

  /**
   * Adds a data node called "name" to the tree, and inserts "data".
   * The data node is added at the current "directory" of this iterator
   * object, so remember to "cd" to the desired location in the tree!
   *
   */
  template <class T>
  PHBoolean AddIODataNode(T* data, const char* name);
};

inline PHDataNodeIterator::PHDataNodeIterator(PHCompositeNode* node)
  : PHNodeIterator(node)
{
}

template <class T>
PHIODataNode<T>*
PHDataNodeIterator::FindIODataNode(PHIODataNode<T>* node,
                                   const char* name)
{
  // TODO:  also check that "name" is not a null string!
  if (!name)
  {
    return 0;
  }
  // Can't do dynamic_cast here; it fails if node was created as
  // PHIODataNode<X> instead of PHIODataNode<T>, even if T is a
  // derived class of X!
  // In general, T -> X does not imply A<T> -> A<X>.
  // ("->" denotes "derives from", and "A" is any template class)
  node = static_cast<PHIODataNode<T>*>(findFirst("PHIODataNode",
                                                 name));
  return node;
}

template <class T>
PHBoolean
PHDataNodeIterator::AddIODataNode(T* data, const char* name)
{
  // TODO:  also check that "name" is not a null string!
  if (!name)
  {
    return false;
  }
  // For IODataNode, ought to check (if possible) that T derives
  // from TObject.  Will typeid() give us this info?

  PHIODataNode<T>* n = new PHIODataNode<T>(data, name);
  if (!n)
  {
    return false;  // problem creating node?
  }
  return addNode(n);
}

#endif
