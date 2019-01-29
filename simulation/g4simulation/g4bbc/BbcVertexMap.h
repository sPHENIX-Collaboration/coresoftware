#ifndef G4BBC_BBCVERTEXMAP_H
#define G4BBC_BBCVERTEXMAP_H

#include "BbcVertex.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class BbcVertexMap : public PHObject {
  
public:

  typedef std::map<unsigned int, BbcVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, BbcVertex*>::iterator            Iter;
  
  virtual ~BbcVertexMap() {}

  virtual void identify(std::ostream &os = std::cout) const {os << "BbcVertexMap base class" << std::endl;}
  virtual void Reset()                                      {}
  virtual int  isValid() const                              {return 0;}
  
  virtual bool   empty()                   const            {return true;}
  virtual size_t size()                    const            {return 0;}
  virtual size_t count(unsigned int idkey) const            {return 0;}
  virtual void   clear()                                    {}
  
  virtual const BbcVertex* get(unsigned int idkey) const    {return NULL;}
  virtual       BbcVertex* get(unsigned int idkey)          {return NULL;}
  virtual       BbcVertex* insert(BbcVertex* vertex)        {return NULL;}
  virtual       size_t     erase(unsigned int idkey)        {return 0;}

  virtual ConstIter begin()                   const {return std::map<unsigned int,BbcVertex*>().end();}
  virtual ConstIter  find(unsigned int idkey) const {return std::map<unsigned int,BbcVertex*>().end();}
  virtual ConstIter   end()                   const {return std::map<unsigned int,BbcVertex*>().end();}

  virtual Iter begin()                   {return std::map<unsigned int,BbcVertex*>().end();}
  virtual Iter  find(unsigned int idkey) {return std::map<unsigned int,BbcVertex*>().end();}
  virtual Iter   end()                   {return std::map<unsigned int,BbcVertex*>().end();}

protected:
  BbcVertexMap() {}

private:

  ClassDef(BbcVertexMap, 1);
};

#endif // G4BBC_BBCVERTEXMAP_H
