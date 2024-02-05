#ifndef TRACKBASEHISTORIC_SVTXTRACKINFO_H
#define TRACKBASEHISTORIC_SVTXTRACKINFO_H

#include <trackbase/TrkrDefs.h>
#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>

#include <limits.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <cstdint>

class SvtxTrackInfo : public PHObject
{
 public:

  ~SvtxTrackInfo() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackInfo base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const SvtxTrackInfo& ) 
  {}

  //! copy content from base class
  virtual void CopyFrom( SvtxTrackInfo* ) 
  {}

  //
  // basic track information ---------------------------------------------------
  //


  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float) {}

  virtual uint8_t get_ndf() const { return uint8_t(-1); }
  virtual void set_ndf(uint8_t) {}

  //bitmap for clusters
  virtual uint64_t get_hitbitmap() const { return uint64_t(-1); }
  virtual void set_hitbitmap(uint64_t) {}


  //vertex information
  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_pos(unsigned int) const { return NAN; }

  virtual float get_px() const { return NAN; }
  virtual void set_px(float) {}

  virtual float get_py() const { return NAN; }
  virtual void set_py(float) {}

  virtual float get_pz() const { return NAN; }
  virtual void set_pz(float) {}

  virtual float get_mom(unsigned int) const { return NAN; }

  virtual float get_p() const { return NAN; }
  virtual float get_pt() const { return NAN; }
  virtual float get_eta() const { return NAN; }
  virtual float get_phi() const { return NAN; }

  virtual float get_covariance(int /*i*/) const { return NAN; }
  virtual void set_covariance(int /*i*/, float /*value*/) {}
  
  virtual short int get_crossing() const { return SHRT_MAX; }
  virtual void set_crossing(short int /*crossing*/) {}
 

private:
  ClassDefOverride(SvtxTrackInfo, 1);
};

#endif