#ifndef __PDBMAPINTVECFLOAT_HH__
#define __PDBMAPINTVECFLOAT_HH__

#include "PdbCalChan.hh"

#include <map>
#include <vector>

class PdbMapIntVecFloat : public PdbCalChan 
{
public:
  PdbMapIntVecFloat() {}
  virtual ~PdbMapIntVecFloat();

  void get_vector(int key, std::vector<float> &fvec) const;
  int add_vector(int key, std::vector<float> &fvec);
  void Verbosity(const int i);
  void Clear(Option_t* ="");
  virtual void print() const;

private:
  std::map<int, std::vector<float> > TheMap;

  ClassDef(PdbMapIntVecFloat,1);

};

#endif /* __PDBMAPINTVECFLOAT_HH__ */
