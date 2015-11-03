#ifndef __PDBMAPINT_HH__
#define __PDBMAPINT_HH__

#include "PdbCalChan.hh"

#include <map>
#include <vector>

class PdbMapIntInt : public PdbCalChan 
{
public:
  PdbMapIntInt(); 
  virtual ~PdbMapIntInt();

  std::map<int,int> get_map() { return TheMap; }
  void set_map(std::map<int,int> *inmap);
  void Verbosity(const int i) { verbosity = i; }
  void Clear(Option_t* ="");
  virtual void print() const;

private:
  std::map<int,int> TheMap;
  int verbosity;

  ClassDef(PdbMapIntInt,1);

};

#endif 
