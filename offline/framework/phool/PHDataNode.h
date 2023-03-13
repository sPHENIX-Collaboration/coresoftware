#ifndef PHOOL_PHDATANODE_H
#define PHOOL_PHDATANODE_H

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
  ~PHDataNode() override;

 public:
  T* getData() { return data.data; }
  void setData(T* d) { data.data = d; }
  void prune() override {}
  void forgetMe(PHNode*) override {}
  void print(const std::string&) override;
  bool write(PHIOManager*, const std::string& = "") override
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
                          const std::string& n)
  : PHNode(n)
{
  type = "PHDataNode";
  setData(d);
}

template <class T>
PHDataNode<T>::PHDataNode(T* d,
                          const std::string& n,
                          const std::string& objtype)
  : PHNode(n, objtype)
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

#endif
