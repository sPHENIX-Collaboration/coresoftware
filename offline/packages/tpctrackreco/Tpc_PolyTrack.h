#pragma once

#include <phool/PHObject.h>

#include <iostream>

class Tpc_PolyTrack : public PHObject
{
 public:
  Tpc_PolyTrack() = default;
  ~Tpc_PolyTrack() override = default;

  void identify(std::ostream& os = std::cout) const override
  {
    os << "Tpc_PolyTrack base class" << std::endl;
  }
  void Reset() override {}
  int isValid() const override { return 0; }

  virtual unsigned int get_event() const { return 0; }
  virtual unsigned int get_track_id() const { return 0; }
  virtual unsigned int get_source_assembled_track_id() const { return 0; }
  virtual int get_fit_status() const { return 0; }
  virtual unsigned int get_nclusters() const { return 0; }
  virtual double get_x() const { return 0.0; }
  virtual double get_y() const { return 0.0; }
  virtual double get_z() const { return 0.0; }
  virtual double get_px() const { return 0.0; }
  virtual double get_py() const { return 0.0; }
  virtual double get_pz() const { return 0.0; }
  virtual double get_charge() const { return 0.0; }
  virtual double get_chi2() const { return 0.0; }
  virtual double get_ndf() const { return 0.0; }
  virtual double get_dedx() const { return 0.0; }
  virtual double get_cov(unsigned int, unsigned int) const { return 0.0; }

  virtual void set_event(unsigned int) {}
  virtual void set_track_id(unsigned int) {}
  virtual void set_source_assembled_track_id(unsigned int) {}
  virtual void set_fit_status(int) {}
  virtual void set_nclusters(unsigned int) {}
  virtual void set_x(double) {}
  virtual void set_y(double) {}
  virtual void set_z(double) {}
  virtual void set_px(double) {}
  virtual void set_py(double) {}
  virtual void set_pz(double) {}
  virtual void set_charge(double) {}
  virtual void set_chi2(double) {}
  virtual void set_ndf(double) {}
  virtual void set_dedx(double) {}
  virtual void set_cov(unsigned int, unsigned int, double) {}

 private:
  ClassDefOverride(Tpc_PolyTrack, 0)
};
