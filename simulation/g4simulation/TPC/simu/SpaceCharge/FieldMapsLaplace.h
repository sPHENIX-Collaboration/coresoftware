#ifndef __FIELDMAPSLAPLACE_H__
#define __FIELDMAPSLAPLACE_H__

//===========================================================
/// \file FieldMapsLaplace.h
/// \brief Implementation of E distortions from Laplace solutions
/// \author Carlos Perez Lara
//===========================================================

#include "FieldMaps.h"

class FieldMapsLaplace:public FieldMaps {
 public:
  virtual void ComputeE();
};

#endif /* __FIELDMAPSLAPLACE_H__ */
