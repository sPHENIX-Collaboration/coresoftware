#ifndef PHOOL_PHIODATANODE_H
#define PHOOL_PHIODATANODE_H

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
class PHIODataNode : public PHDataNode<T>
{
  friend class PHNodeIOManager;

 public:
  T *operator*() { return this->getData(); }
  PHIODataNode(T *, const std::string &);
  PHIODataNode(T *, const std::string &, const std::string &);
  virtual ~PHIODataNode() {}
  typedef PHTypedNodeIterator<T> iterator;
  void BufferSize(int size) {buffersize = size;}
  void SplitLevel(int split) {splitlevel = split;}

 protected:
  bool write(PHIOManager *, const std::string & = "") override;
  PHIODataNode() = delete;
  int buffersize;
  int splitlevel;
};

template <class T>
PHIODataNode<T>::PHIODataNode(T *d, const std::string &n)
  : PHDataNode<T>(d, n)
  , buffersize(32000)
  , splitlevel(99)
{
  this->type = "PHIODataNode";
  TObject *TO = static_cast<TObject *>(d);
  this->objectclass = TO->GetName();
}

template <class T>
PHIODataNode<T>::PHIODataNode(T *d, const std::string &n,
                              const std::string &objtype)
  : PHDataNode<T>(d, n, objtype)
  , buffersize(32000)
  , splitlevel(99)
{
  this->type = "PHIODataNode";
  TObject *TO = static_cast<TObject *>(d);
  this->objectclass = TO->GetName();
}

template <class T>
bool PHIODataNode<T>::write(PHIOManager *IOManager, const std::string &path)
{
  if (this->persistent)
  {
    PHNodeIOManager *np = dynamic_cast<PHNodeIOManager *>(IOManager);
    if (np)
    {
      std::string newPath = path + phooldefs::branchpathdelim + this->name;
      bool bret = false;
      if (dynamic_cast<TObject *>(this->data.data))
      {
        bret = np->write(&(this->data.tobj), newPath, buffersize, splitlevel);
      }
      return bret;
    }
  }
  return true;
}

#endif /* PHOOL_PHIODATANODE_H */
