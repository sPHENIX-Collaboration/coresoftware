// $Id: PHG4Hitv6.h,v 1.2 2015/01/06 02:52:08 jinhuang Exp $                                                                                             

/*!
 * \file PHG4Hitv6.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.2 $
 * \date $Date: 2015/01/06 02:52:08 $
 */

#ifndef PHG4HITV6_H_
#define PHG4HITV6_H_

#include "PHG4Hitv5.h"
#include <Geant4/G4Allocator.hh>

#ifndef __CINT__
class PHG4Hitv6;
extern G4Allocator<PHG4Hitv6> PHG4Hitv6Allocator;
#endif

/*!
 * \brief PHG4Hitv6 with scintillation light yield and pathlength
 */
class PHG4Hitv6 : public PHG4Hitv5
{
public:
  PHG4Hitv6();
  PHG4Hitv6(const PHG4Hit& g4hit);
  virtual
  ~PHG4Hitv6();
  virtual void
  print() const;

  float get_light_yield() const
    {
      return light_yield;
    }

  void
  set_light_yield(float lightYield)
  {
    light_yield = lightYield;
  }

  float
  get_path_length() const
  {
    return path_length;
  }

  void
  set_path_length(float pathLength)
  {
    path_length = pathLength;
  }

protected:
  //! a number proportional to the scintillation light yield.
  float light_yield;

  //! path length of the track to the hit
  float path_length;

  ClassDef(PHG4Hitv6,1)

};

#endif /* PHG4HITV6_H_ */
