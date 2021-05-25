// Tell emacs that this is a C++ source
//  -*- C++ -*-.
// $Id: $

/*!
 * \file PHG4HitEval.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef G4MAIN_PHG4HITEVAL_H
#define G4MAIN_PHG4HITEVAL_H

#include "PHG4Hitv1.h"

class PHG4Hit;
class PHObject;

/*!
 * \brief PHG4HitEval for evaluating PHG4Hitv1 in CINT readable evaluation trees
 */
class PHG4HitEval : public PHG4Hitv1
{
 public:
  PHG4HitEval();

  PHG4HitEval(const PHG4Hit *g4hit);

  ~PHG4HitEval() override {}

  void CopyFrom(const PHObject *phobj) override;

  float
  get_eion() const override
  {
    return eion;
  }
  void
  set_eion(const float f) override
  {
    eion = f;
  }

  int get_scint_id() const override
  {
    return scint_id;
  }

  void
  set_scint_id(const int i) override
  {
    scint_id = i;
  }

  float
  get_light_yield() const override
  {
    return light_yield;
  }

  void
  set_light_yield(float lightYield) override
  {
    light_yield = lightYield;
  }

  float
  get_path_length() const override
  {
    return path_length;
  }

  void
  set_path_length(float pathLength) override
  {
    path_length = pathLength;
  }

 protected:
  float eion;

  int scint_id;

  //! a number proportional to the scintillation light yield.
  float light_yield;

  //! path length of the track to the hit
  float path_length;

  ClassDefOverride(PHG4HitEval, 1)
};

#endif /* PHG4HITEVAL_H_ */
