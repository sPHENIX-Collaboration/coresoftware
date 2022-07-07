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

#include <TGeoMatrix.h>

#include <cmath>
#include <iostream>

class ActsGeometry;
class TVector2;
class TVector3;
class PHG4Hit;

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

  //! print information about this layer
  void identify(std::ostream&) const override;

  //!@name accessors
  //@{
  int get_layer() const override {return m_layer;}
  double get_radius() const override {return m_radius;}
  double get_thickness() const override { return m_thickness;}
  double get_zmin() const override {return m_zmin;}
  double get_zmax() const override {return m_zmax;}
  double get_pitch() const { return m_pitch; }

  //! segmentation type
  MicromegasDefs::SegmentationType get_segmentation_type() const {return m_segmentation_type;}

  //! drift direction
  MicromegasDefs::DriftDirection get_drift_direction() const {return m_drift_direction;}

  // check if hit radius matches this cylinder
  bool check_radius( const TVector3& ) const;

  //! convert world to local position coordinates in (planar) tile reference frame
  TVector3 get_local_from_world_coords( uint tileid, ActsGeometry*, const TVector3& ) const;

  //! convert world to local direction coordinates in (planar) tile reference frame
  TVector3 get_local_from_world_vect( uint tileid, ActsGeometry*, const TVector3& ) const;

  //! convert local to world position coordinates in (planar) tile reference frame
  TVector3 get_world_from_local_coords( uint tileid, ActsGeometry*, const TVector2& ) const;

  //! convert local to world position coordinates in (planar) tile reference frame
  TVector3 get_world_from_local_coords( uint tileid, ActsGeometry*, const TVector3& ) const;

  //! convert local to world direction coordinates in (planar) tile reference frame
  TVector3 get_world_from_local_vect( uint tileid, ActsGeometry*, const TVector3& ) const;

  //! get tile for a given world location assuming tiles are portion of cylinder centered around tile center
  /** it is used to find the tile hit by a given track.
   * The track is first extrapolated to the layer cylinder, the relevant tile is found, if any
   * the track is then extrapolated a second time to the correct tile plane
   */
  int find_tile_cylindrical( const TVector3& ) const;

  //! get number of tiles
  size_t get_tiles_count() const { return m_tiles.size(); }

  //! get tile for given tileid
  const MicromegasTile& get_tile( uint tileid ) const
  {
    assert( tileid < m_tiles.size() );
    return m_tiles[tileid];
  }

  //! get strip for a give world location and tile
  int find_strip_from_world_coords( uint tileid, ActsGeometry*, const TVector3& ) const;

  //! get strip for a give world location and tile
  int find_strip_from_local_coords( uint tileid, ActsGeometry*, const TVector2& ) const;

  //! get strip length for a given tile
  double get_strip_length( uint tileid, ActsGeometry* ) const;

  //! get number of strips
  uint get_strip_count( uint tileid, ActsGeometry* ) const;

  //! get local coordinates for a given tile and strip
  TVector2 get_local_coordinates( uint tileid, ActsGeometry*, uint stripnum ) const;

  //! get world coordinates for a given tile and strip
  TVector3 get_world_coordinates( uint tileid, ActsGeometry*, uint stripnum ) const;

  /// reference radius used in macro to convert tile size in azimuth into a angle range (cm)
  static constexpr double reference_radius = 82;

  //@}

  //!@name modifiers
  //@{
  void set_layer(const int i) override {m_layer = i;}
  void set_radius(const double value) override {m_radius = value;}
  void set_thickness(const double value) override {m_thickness = value;}
  void set_zmin(const double value) override {m_zmin = value;}
  void set_zmax(const double value) override {m_zmax = value;}
  void set_pitch( double value ) { m_pitch = value; }

  //! tiles
  void set_tiles( const MicromegasTile::List& tiles ) {m_tiles = tiles;}

  //! segmentation
  void set_segmentation_type( MicromegasDefs::SegmentationType value ) {m_segmentation_type = value;}

  //! drift direction
  void set_drift_direction( MicromegasDefs::DriftDirection value ) {m_drift_direction = value;}
  //@}

  private:

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
  /** 
   * \brief tiles are only used in "find_tile_cylindrical". 
   * For all other methods we use ACTS surfaces instead 
   */
  MicromegasTile::List m_tiles;

  ClassDefOverride(CylinderGeomMicromegas, 1)
};

#endif
