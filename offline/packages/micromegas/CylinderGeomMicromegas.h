// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef MICROMEGAS_CYLINDERGEOMMICROMEGAS_H
#define MICROMEGAS_CYLINDERGEOMMICROMEGAS_H

/*!
 * \file CylinderGeomMicromegas.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasDefs.h"
#include "MicromegasTile.h"
#include <g4detectors/PHG4CylinderGeom.h>

#include <cmath>
#include <iostream>

class TVector3;

class CylinderGeomMicromegas : public PHG4CylinderGeom
{
  public:

  //* empty constructor
  explicit CylinderGeomMicromegas()
  {}

  //* constructor
  CylinderGeomMicromegas(int layer, MicromegasDefs::SegmentationType type )
    : m_layer( layer )
    , m_segmentation_type( type )
  {}

  //!@name accessors
  //@{
  virtual int get_layer() const {return m_layer;}
  virtual double get_radius() const {return m_radius;}
  virtual double get_thickness() const { return m_thickness;}
  virtual double get_zmin() const {return m_zmin;}
  virtual double get_zmax() const {return m_zmax;}
  double get_pitch() const { return m_pitch; }

  //! segmentation type
  MicromegasDefs::SegmentationType get_segmentation_type() const {return m_segmentation_type;}

  //! get tile and strip for a give world location
  std::pair<int,int> find_strip( const TVector3& ) const;

  //! get strip length for a given tile
  double get_strip_length( uint tileid ) const;

  //! get world location for a given tile and strip
  TVector3 get_world_coordinate( uint tileid, uint stripnum ) const;

  //! print information about this layer
  virtual void identify(std::ostream&) const;

  //@}

  //!@name modifiers
  //@{
  virtual void set_layer(const int i) {m_layer = i;}
  virtual void set_radius(const double value) {m_radius = value;}
  virtual void set_thickness(const double value) {m_thickness = value;}
  virtual void set_zmin(const double value) {m_zmin = value;}
  virtual void set_zmax(const double value) {m_zmax = value;}
  void set_pitch( double value ) { m_pitch = value; }

  //! tiles
  void set_tiles( const MicromegasTile::List& tiles ) { m_tiles = tiles; }
  void set_segmentation_type( MicromegasDefs::SegmentationType value ) {m_segmentation_type = value;}
  //@}

  private:

  int m_layer = 0;
  MicromegasDefs::SegmentationType m_segmentation_type = MicromegasDefs::SegmentationType::SEGMENTATION_PHI;
  double m_radius = 0;
  double m_thickness = 0;
  double m_zmin = 0;
  double m_zmax = 0;

  // 1mm pitch by default
  double m_pitch = 0.1;

  //! tiles
  MicromegasTile::List m_tiles;

  ClassDef(CylinderGeomMicromegas, 1)
};

#endif
