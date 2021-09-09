#ifndef __SVTXTRACKSTATE_H__
#define __SVTXTRACKSTATE_H__

#include <phool/PHObject.h>
#include <cmath>

class SvtxTrackState : public PHObject
{
 public:
  ~SvtxTrackState() override {}
    
  void identify(std::ostream &os = std::cout) const override
  {
    os << "SvtxTrackState base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject *CloneMe() const override { return nullptr; }

  virtual float get_pathlength() const { return NAN; }

  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_pos(unsigned int /*i*/) const { return NAN; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float) {}

  virtual float get_py() const { return NAN; }
  virtual void set_py(float) {}

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float) {}

  virtual float get_mom(unsigned int /*i*/) const { return NAN; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }

  virtual float get_error(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void set_error(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}

  virtual std::string get_name() const { return ""; }
  virtual void set_name(const std::string &/*name*/) {}

  ///@name convenience interface, also found in Trkrcluster
  //@{

  /// rphi error
  virtual float get_rphi_error() const
  {
    return NAN;
  }

  /// phi error
  virtual float get_phi_error() const
  {
    return NAN;
  }

  /// z error
  virtual float get_z_error() const
  {
    return NAN;
  }

  //@}

 protected:
  SvtxTrackState(float /*pathlength*/ = 0.0) {}

  ClassDefOverride(SvtxTrackState, 1);
};

#endif
