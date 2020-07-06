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
  CylinderGeomMicromegas(int layer )
    : m_layer( layer )
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

  //! drift direction
  MicromegasDefs::DriftDirection get_drift_direction() const {return m_drift_direction;}

  //! get tile for a given world location
  int find_tile( const TVector3& ) const;

  //! get tile and strip for a give world location
  std::pair<int,int> find_strip( const TVector3& ) const;

  //! get strip for a give world location and tile
  int find_strip( uint tileid, const TVector3& ) const;

  //! get strip length for a given tile
  double get_strip_length( uint tileid ) const;

  //! get number of strips
  uint get_strip_count( uint tileid ) const;

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
  void set_tiles( const MicromegasTile::List& tiles ) {m_tiles = tiles;}

  //! segmentation
  void set_segmentation_type( MicromegasDefs::SegmentationType value ) {m_segmentation_type = value;}

  //! drift direction
  void set_drift_direction( MicromegasDefs::DriftDirection value ) {m_drift_direction = value;}
  //@}

  private:

  // check if hit radius matches this cylinder
  bool check_radius( const TVector3& ) const;

  //! layer id
  int m_layer = 0;

  //! segmentation type
  MicromegasDefs::SegmentationType m_segmentation_type = MicromegasDefs::SegmentationType::SEGMENTATION_PHI;

  //! drift direction
  MicromegasDefs::DriftDirection m_drift_direction = MicromegasDefs::DriftDirection::OUTWARD;

  //! layer radius
  double m_radius = 0;

  //! layer thickness
  double m_thickness = 0;

  //! layer z extend
  double m_zmin = 0;

  //! layer z extend
  double m_zmax = 0;

  //! 1mm pitch by default
  double m_pitch = 0.1;

  //! tiles
  MicromegasTile::List m_tiles;

  ClassDef(CylinderGeomMicromegas, 1)
};

#endif
