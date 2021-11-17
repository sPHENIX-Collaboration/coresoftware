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
  /**
   * each (planar) tile has a local ref system defined as:
   * - origin at center of the tile
   * - z axis same as phenix,
   * - y axis perpendicular to the surface, outward,
   * - x axis perpendicular to y and z to have a direct ref. frame
   **/
  TVector3 get_local_from_world_coords( uint tileid, const TVector3& ) const;

  //! convert world to local direction coordinates in (planar) tile reference frame
  TVector3 get_local_from_world_vect( uint tileid, const TVector3& ) const;

  //! convert local to world position coordinates in (planar) tile reference frame
  /**
   * each (planar) tile has a local ref system defined as:
   * - origin at center of the tile
   * - z axis same as phenix,
   * - y axis perpendicular to the surface, outward,
   * - x axis perpendicular to y and z to have a direct ref. system
   **/
  TVector3 get_world_from_local_coords( uint tileid, const TVector3& ) const;

  //! convert local to world direction coordinates in (planar) tile reference frame
  TVector3 get_world_from_local_vect( uint tileid, const TVector3& ) const;

  //! get tile for a given world location assuming tiles are portion of cylinder centered around tile center
  int find_tile_cylindrical( const TVector3& ) const;

  //! get tile for a given world location assuming tiles are planes centered on tile center and tengent to local cylinder
  int find_tile_planar( const TVector3& ) const;

  //! get number of tiles
  size_t get_tiles_count() const { return m_tiles.size(); }

  //! get tile for given tileid
  const MicromegasTile& get_tile( uint tileid ) const
  {
    assert( tileid < m_tiles.size() );
    return m_tiles[tileid];
  }

  //! convert g4hit coordinates from cylinder Micromegas to planar
  /**
    * this assumes that Micromegas Geant4 implementation are cylinders, while actual tiles are planes
    * one must then 'drift' the g4hit along its momentum from its radius to the releval "y" in local coordinates
    */
  void convert_to_planar( uint tileid, PHG4Hit* ) const;

  //! get strip for a give world location and tile
  int find_strip_from_world_coords( uint tileid, const TVector3& ) const;

  //! get strip for a give world location and tile
  int find_strip_from_local_coords( uint tileid, const TVector3& ) const;

  //! get strip length for a given tile
  double get_strip_length( uint tileid ) const;

  //! get number of strips
  uint get_strip_count( uint tileid ) const;

  //! get local coordinates for a given tile and strip
  TVector3 get_local_coordinates( uint tileid, uint stripnum ) const;

  //! get world coordinates for a given tile and strip
  TVector3 get_world_coordinates( uint tileid, uint stripnum ) const;

  //! get phi angle at center of tile
  double get_center_phi( uint tileid ) const
  {
    assert( tileid < m_tiles.size() );
    return m_tiles[tileid].m_centerPhi;
  }

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

  // get local to master transformation matrix for a given tile
  TGeoHMatrix transformation_matrix( uint tileid ) const;

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

  ClassDefOverride(CylinderGeomMicromegas, 1)
};

#endif
