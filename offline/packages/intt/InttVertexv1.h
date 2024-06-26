#ifndef INTT_INTTVERTEXV1_H
#define INTT_INTTVERTEXV1_H

#include "InttVertex.h"

#include <iostream>

class PHObject;

class InttVertexv1 : public InttVertex
{
 public:
  InttVertexv1();
  ~InttVertexv1() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = InttVertexv1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new InttVertexv1(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_x() const override { return _pos[0]; }
  void set_x(float x) override { _pos[0] = x; }

  float get_y() const override { return _pos[1]; }
  void set_y(float y) override { _pos[1] = y; }

  float get_z() const override { return _pos[2]; }
  void set_z(float z) override { _pos[2] = z; }

  float get_position(unsigned int coor) const override { return _pos[coor]; }
  void set_position(unsigned int coor, float xi) override { _pos[coor] = xi; }

  float get_error(unsigned int i, unsigned int j) const override;        //< get vertex error covar
  void set_error(unsigned int i, unsigned int j, float value) override;  //< set vertex error covar

  double get_chi2ndf() const override { return chi2ndf; }
  void set_chi2ndf(double val) override { chi2ndf = val; }

  double get_width() const override { return width; }
  void set_width(double val) override { width = val; }

  bool get_good() const override { return good; }
  void set_good(bool val) override { good = val; }

  unsigned int get_nclus() const override { return nclus; }
  void set_nclus(unsigned int val) override { nclus = val; }

  unsigned int get_ntracklet() const override { return ntracklet; }
  void set_ntracklet(unsigned int val) override { ntracklet = val; }

  unsigned int get_ngroup() const override { return ngroup; }
  void set_ngroup(unsigned int val) override { ngroup = val; }

  double get_peakratio() const override { return peakratio; }
  void set_peakratio(double val) override { peakratio = val; }

  double get_peakwidth() const override { return peakwidth; }
  void set_peakwidth(double val) override { peakwidth = val; }

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

 private:
  unsigned int _id = std::numeric_limits<unsigned int>::max();  //< unique identifier within container
  float _pos[3] = {};                                           //< collision position x,y,z
  float _err[6] = {};                                           //< error covariance matrix (+/- cm^2)

  double chi2ndf{-1};
  double width{-1};
  bool good{false};
  unsigned int nclus{0};
  unsigned int ntracklet{0};
  unsigned int ngroup{0};
  double peakratio{0};
  double peakwidth{0};

  ClassDefOverride(InttVertexv1, 1);
};

#endif
