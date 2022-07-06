// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4MICROMEGASDETECTOR_H
#define PHG4MICROMEGASDETECTOR_H

/*!
 * \file PHG4MicromegasDetector.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4main/PHG4Detector.h>
#include <micromegas/MicromegasTile.h>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4MicromegasDisplayAction;
class PHG4Subsystem;
class PHParameters;

class PHG4MicromegasDetector : public PHG4Detector
{
 public:

  //! constructor
  PHG4MicromegasDetector(PHG4Subsystem*, PHCompositeNode*, PHParameters*, const std::string&);

  //! construct
  void ConstructMe(G4LogicalVolume*) override;

  void Print(const std::string &what = "ALL") const override;

  //! set first layer number
  void set_first_layer( int layer ) { m_FirstLayer = layer; }

  //! get first layer number
  int get_first_layer() const { return m_FirstLayer; }

  //! returns true if passed volume is an active volume of this detector
  bool IsInDetector(G4VPhysicalVolume*) const;

  //! return layer associated to a given volume, or -1 if invalid
  int get_layer(G4VPhysicalVolume*) const;

  //! return tile id associated to a given volume, or -1 if invalid
  int get_tileid(G4VPhysicalVolume*) const;

  //! super detector name
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }

  //! super detector name
  const std::string SuperDetector() const { return m_SuperDetector; }

  //! access the display action
  PHG4MicromegasDisplayAction* GetDisplayAction() { return m_DisplayAction; }

  private:

  //! setup tiles
  /** the method is now private because tiles are now hard coded */
  void setup_tiles();
  
  //! create needed material
  void create_materials() const;

  //! construct
  void construct_micromegas(G4LogicalVolume*);

  //! create a micromegas tile of given dimension
  /** returns the master logical volume that can then be placed inside the world logical volume */
  G4LogicalVolume* construct_micromegas_tile( int tileid );
  
  //! construct FEE board
  G4LogicalVolume* construct_fee_board( int id );
  
  //! add geometry node
  /*! this handles the internal (module/strips) segmentation, needed for tracking*/
  void add_geometry_node();
  
  //! vis attribute handling (save memory in batch)
  PHG4MicromegasDisplayAction* m_DisplayAction = nullptr;

  //! detector parameters
  PHParameters* m_Params = nullptr;

  //! map layer index to radius (cm)
  /** it is filled while creating G4 volumes */
  std::map<int, double> m_layer_radius;
  
  //! map layer index to thickness (cm)
  /** it is filled while creating G4 volumes */
  std::map<int, double> m_layer_thickness;

  //! active volumes, and mapping to layer
  /*! it is needed in the stepping action to map a volume to a given layer */
  std::map<G4VPhysicalVolume*, int> m_activeVolumes;

  //! map active volumes to tile number
  /*! it is needed in the stepping action to map a volume to a given tile id */
  std::map<G4VPhysicalVolume*, int> m_tiles_map;

  //! also store passive volumes
  std::set<G4VPhysicalVolume*> m_passiveVolumes;
  
  //! super detector name
  std::string m_SuperDetector;

  //! micromegas tiles
  MicromegasTile::List m_tiles;

  //! first layer number
  /* there are two layers in the detector */
  int m_FirstLayer = 0;

};

#endif // MICROMEGASDETECTOR_H
