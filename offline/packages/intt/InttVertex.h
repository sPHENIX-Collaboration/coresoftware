#ifndef INTT_INTTVERTEX_H
#define INTT_INTTVERTEX_H

#include <globalvertex/Vertex.h>

#include <cmath>
#include <iostream>

class InttVertex : public Vertex
{
 public:
  ~InttVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "InttVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const override { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) override {}

  virtual float get_x() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_x(float) override {}

  virtual float get_y() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_y(float) override {}

  virtual float get_z() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_z(float) override {}

  virtual float get_position(unsigned int) const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_position(unsigned int, float) override {}

  virtual float get_error(unsigned int, unsigned int) const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_error(unsigned int, unsigned int, float) override {}

  virtual double get_chi2ndf() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual void set_chi2ndf(double) {}

  virtual double get_width() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual void set_width(double) {}

  virtual bool get_good() const { return false; }
  virtual void set_good(bool) {}

  virtual unsigned int get_nclus() const { return 0; }
  virtual void set_nclus(unsigned int) {}

  virtual unsigned int get_ntracklet() const { return 0; }
  virtual void set_ntracklet(unsigned int) {}

  virtual unsigned int get_ngroup() const { return 0; }
  virtual void set_ngroup(unsigned int) {}

  virtual double get_peakratio() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual void set_peakratio(double) {}

  virtual double get_peakwidth() const { return std::numeric_limits<double>::quiet_NaN(); }
  virtual void set_peakwidth(double) {}

 protected:
  InttVertex() {}

 private:
  ClassDefOverride(InttVertex, 1);
};

#endif
