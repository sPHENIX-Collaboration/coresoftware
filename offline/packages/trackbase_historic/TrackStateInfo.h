#ifndef TRACKBASEHISTORIC_TRACKSTATEINFO_H
#define TRACKBASEHISTORIC_TRACKSTATEINFO_H

//#include "SvtxTrackState.h"
//#include "TrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <g4main/PHG4HitDefs.h>
#include <phool/PHObject.h>

#include <limits.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

class TrackStateInfo : public PHObject
{
 public:

  TrackStateInfo() = default;
  ~TrackStateInfo() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackStateInfo base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  
  //! copy content from base class
  virtual void CopyFrom( const TrackStateInfo& ) 
  {}

  //! copy content from base class
  virtual void CopyFrom( TrackStateInfo* ) 
  {}

  //
  // basic track information ---------------------------------------------------
  //

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

  //virtual float get_error(int /*i*/, int /*j*/) const { return NAN; }
  //virtual void set_error(int /*i*/, int /*j*/, float /*value*/) {}
  
  virtual float get_covariance(int) const {return NAN;}
  virtual void set_covariance(int, float) {}

 

 private:
  ClassDefOverride(TrackStateInfo, 1);
};

#endif