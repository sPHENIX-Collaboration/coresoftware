#ifndef INTT_INTTVERTEXMAP_H
#define INTT_INTTVERTEXMAP_H

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class InttVertex;

class InttVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, InttVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, InttVertex*>::iterator Iter;

  ~InttVertexMap() override {}

  void identify(std::ostream& os = std::cout) const override { os << "InttVertexMap base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const InttVertex* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual InttVertex* get(unsigned int /*idkey*/) { return nullptr; }
  virtual InttVertex* insert(InttVertex* /*vertex*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  InttVertexMap() {}

 private:
  ClassDefOverride(InttVertexMap, 1);
};

#endif  // INTT_INTTVERTEXMAP_H
