#ifndef GETCLASS_H__
#define GETCLASS_H__

#include "PHDataNode.h"
#include "PHIODataNode.h"
#include "PHNode.h"
#include "PHNodeIterator.h"

#include <TObject.h>

#include <string>

class PHCompositeNode;

namespace findNode
{
template <class T>
T *getClass(PHCompositeNode *top, const std::string &name)
{
  PHNodeIterator iter(top);
  PHNode *FoundNode = iter.findFirst(name.c_str());  // returns pointer to PHNode
  if (!FoundNode)
  {
    return NULL;
  }
  // first test if it is a PHDataNode
  PHDataNode<T> *DNode = dynamic_cast<PHDataNode<T> *>(FoundNode);
  if (DNode)
  {
    T *object = dynamic_cast<T *>(DNode->getData());
    if (object)
    {
      return object;
    }
  }

  // We make a static cast for PHIODataNode<TObject> since all
  // PHIODataNodes have to contain a TObject (otherwise it cannot be
  // written out and it should be a PHDataNode. dynamic cast does not
  // work (see some explanation in PHTypedNodeIterator.h)

  PHIODataNode<TObject> *IONode = static_cast<PHIODataNode<TObject> *>(FoundNode);
  if (IONode)
  {
    T *object = dynamic_cast<T *>(IONode->getData());
    if (!object)
    {
      return NULL;
    }
    else
    {
      return object;
    }
  }

  return NULL;
}
}

#endif /* GETCLASS_H */
