#ifndef __PHIODATANODE_H__
#define __PHIODATANODE_H__

//  Declaration of class PHIODataNode which can hold persistent data
//  Author: Matthias Messer

#include "PHDataNode.h"
#include "PHIOManager.h"
#include "PHNodeIOManager.h"
#include "PHTypedNodeIterator.h"

#include <TObject.h>

class TObject;

template <typename T>
class PHIODataNode : public PHDataNode <T>
{
  friend class PHNodeIOManager;

 public:
  T* operator*() {return this->getData();}
  PHIODataNode(T*, const PHString&);
  PHIODataNode(T*, const PHString&, const PHString&);
  virtual ~PHIODataNode() {}
  typedef PHTypedNodeIterator<T> iterator;

 protected:
  virtual PHBoolean write(PHIOManager *, const PHString& = "");
  PHIODataNode() {}
};

template <class T>
PHIODataNode<T>::PHIODataNode(T* d, const PHString& name)
  : PHDataNode<T>(d, name)
{
  this->type = "PHIODataNode";
}

template <class T>
PHIODataNode<T>::PHIODataNode(T* d, const PHString& name,
                              const PHString& objtype)
  : PHDataNode<T>(d, name, objtype)
{
  this->type = "PHIODataNode";
}

template <class T>
PHBoolean
PHIODataNode<T>::write(PHIOManager* IOManager, const PHString& path)
{
  if (this->persistent)
    {
      PHNodeIOManager *np = dynamic_cast<PHNodeIOManager*>(IOManager);
      if (np)
        {
          PHString newPath = path + "/" + this->name;
	  PHBoolean bret = False;
	  if (dynamic_cast<TObject *> (this->data.data))
	    {
	      bret =  np->write(&(this->data.tobj), newPath);
	    }
	  return bret;
        }
    }
  return True;
}

#endif /* __PHIODATANODE_H__ */


