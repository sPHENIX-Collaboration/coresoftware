#ifndef __PHIODATANODE_H__
#define __PHIODATANODE_H__

//  Declaration of class PHIODataNode which can hold persistent data
//  Author: Matthias Messer

#include "PHDataNode.h"
#include "PHIOManager.h"
#include "PHNodeIOManager.h"
#include "PHTypedNodeIterator.h"
#include "phooldefs.h"

#include <TObject.h>

#include <string>

template <typename T>
class PHIODataNode : public PHDataNode <T>
{
  friend class PHNodeIOManager;

 public:
  T* operator*() {return this->getData();}
  PHIODataNode(T*, const std::string &);
  PHIODataNode(T*, const std::string &, const std::string &);
  virtual ~PHIODataNode() {}
  typedef PHTypedNodeIterator<T> iterator;

 protected:
  virtual PHBoolean write(PHIOManager *, const PHString& = "");
  PHIODataNode() {}
};

template <class T>
PHIODataNode<T>::PHIODataNode(T* d, const std::string& name)
  : PHDataNode<T>(d, name)
{
  this->type = "PHIODataNode";
}

template <class T>
PHIODataNode<T>::PHIODataNode(T* d, const std::string& name,
                              const std::string& objtype)
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
	  std::string tmpstr = path.getString() + phooldefs::branchpathdelim + this->name;
          PHString newPath = tmpstr.c_str();
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
