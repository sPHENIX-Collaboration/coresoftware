#ifndef CYLINDERGEOMMICROMEGAS_H
#define CYLINDERGEOMMICROMEGAS_H

/*!
 * \file CylinderGeomMicromegas.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"

#include <g4detectors/PHG4CylinderGeom.h>

#include <cmath>
#include <iostream>

class CylinderGeomMicromegas : public PHG4CylinderGeom
{
  public:

  //* constructor
  CylinderGeomMicromegas(int layer, MicromegasDefs::SegmentationType type )
    : m_layer( layer )
    , m_segmentation_type( type )
  {}

  //!@name accessors
  //@{
  int get_layer() const {return m_layer;}
  double get_radius() const {return m_radius;}
  double get_thickness() const { return m_thickness;}
  double get_zmin() const {return m_zmin;}
  double get_zmax() const {return m_zmax;}

  MicromegasDefs::SegmentationType get_segmentation_type() const {return m_segmentation_type;}

  //@}

  //!@name reimplemented modifiers
  //@{
  void set_layer(const int i) {m_layer = i;}
  void set_radius(const double value) {m_radius = value;}
  void set_thickness(const double value) {m_thickness = value;}
  void set_zmin(const double value) {m_zmin = value;}
  void set_zmax(const double value) {m_zmax = value;}
  void set_segmentation_type( MicromegasDefs::SegmentationType value ) {m_segmentation_type = value;}
  //@}

  private:

  int m_layer = 0;
  MicromegasDefs::SegmentationType m_segmentation_type = MicromegasDefs::SegmentationType::SEGMENTATION_PHI;
  double m_radius = 0;
  double m_thickness = 0;
  double m_zmin = 0;
  double m_zmax = 0;

  ClassDef(CylinderGeomMicromegas, 1)
};

#endif
