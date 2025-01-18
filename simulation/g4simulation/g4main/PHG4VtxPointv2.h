// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4VTXPOINTV2_H
#define G4MAIN_PHG4VTXPOINTV2_H

#include "PHG4VtxPointv1.h"

#include "PHG4MCProcessDefs.h"

#include <climits>   // for INT_MIN
#include <cmath>     // def of NAN
#include <iostream>  // for cout, ostream

class PHG4VtxPointv2 : public PHG4VtxPointv1
{
 public:
  PHG4VtxPointv2() = default;
  PHG4VtxPointv2(const PHG4VtxPointv2& rhs) = default;
  PHG4VtxPointv2& operator=(const PHG4VtxPointv2&) = default;
  PHG4VtxPointv2(PHG4VtxPointv2&& rhs) = default;
  PHG4VtxPointv2& operator=(PHG4VtxPointv2&&) = default;
  explicit PHG4VtxPointv2(const PHG4VtxPoint* vtx)
    : PHG4VtxPointv1(vtx)
  {
  }

  PHG4VtxPointv2(const double x, const double y, const double z, const double t,
                 const int id_value = std::numeric_limits<int>::min(),
                 const PHG4MCProcess process = PHG4MCProcess::kPNoProcess)
    : PHG4VtxPointv1(x, y, z, t, id_value)
  {
    setProcess(process);
  };

  ~PHG4VtxPointv2() override = default;

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  /// set process property
  void setProcess(int proc)
  {
    auto prop = ((PropEncoding) mProp);
    prop.properties.process = proc;
    mProp = prop.i;
  }

  /// get the production process (id) of this track
  int getProcess() const { return ((PropEncoding) mProp).properties.process; }
  std::string_view getProdProcessAsString() const;

 protected:
  // internal structure to allow convenient manipulation
  // of properties as bits on an int
  union PropEncoding
  {
    explicit PropEncoding(int a)
      : i(a)
    {
    }
    int i;
    struct
    {
      int storage : 1;           // encoding whether to store this track to the output
      unsigned int process : 6;  // encoding process that created this track (enough to store TMCProcess from ROOT)
    } properties;
  };

 private:
  int mProp = 0;

  ClassDefOverride(PHG4VtxPointv2, 1)
};

inline std::string_view PHG4VtxPointv2::getProdProcessAsString() const
{
  auto procID = getProcess();
  if (procID >= 0)
  {
    return PHG4MCProcessName[procID];
  }
  else
  {
    return PHG4MCProcessName[PHG4MCProcess::kPNoProcess];
  }
}

#endif
