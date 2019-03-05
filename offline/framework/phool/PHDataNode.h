#ifndef PHDATANODE_H__
#define PHDATANODE_H__

//  Declaration of class PHDataNode
//  Purpose: a node which can hold a data object (template)

#include "PHNode.h"

#include <iostream>
class TObject;

template <class T>
class PHDataNode : public PHNode
{
 public:
  PHDataNode(T*, const std::string&);
  PHDataNode(T*, const std::string&, const std::string&);
  virtual ~PHDataNode();

 public:
  T* getData() { return data.data; }
  void setData(T* d) { data.data = d; }
  virtual void prune() {}
  virtual void forgetMe(PHNode*) {}
  virtual void print(const std::string&);
  virtual bool write(PHIOManager*, const std::string& = "")
  {
    return true;
  }

 protected:
  union tobjcast {
    T* data;
    TObject* tobj;
  };
  tobjcast data;
  PHDataNode() = delete;
};


template <class T>
PHDataNode<T>::PHDataNode(T* d,
                          const std::string& name)
  : PHNode(name)
{
  type = "PHDataNode";
  setData(d);
}

template <class T>
PHDataNode<T>::PHDataNode(T* d,
                          const std::string& name,
                          const std::string& objtype)
  : PHNode(name, objtype)
{
  type = "PHDataNode";
  setData(d);
}

template <class T>
PHDataNode<T>::~PHDataNode()
{
  // This means that the node has complete responsibility for the
  // data it contains. Check for null pointer just in case some
  // joker adds a node with a null pointer
  if (data.data)
  {
    delete data.data;
    data.data = 0;
  }
}

template <class T>
void PHDataNode<T>::print(const std::string& path)
{
  std::cout << path << name << " (";
  if (!this->objectclass.empty())
  {
    if (type.find("IO") != std::string::npos)
    {
      std::cout << "IO";
    }
    else
    {
      std::cout << "Data";
    }
    std::cout << "," << this->objectclass;
  }
  else
  {
    std::cout << type;
  }
  std::cout << ")" << std::endl;
}

#endif /* __PHDATANODE_H__ */
