// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4MICROMEGASDETECTOR_H
#define PHG4MICROMEGASDETECTOR_H

/*!
 * \file PHG4MicromegasDetector.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
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

  //! super detector name
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }

  //! super detector name
  const std::string SuperDetector() const { return m_SuperDetector; }

  private:

  //! create needed material
  void create_materials() const;

  //! construct
  void construct_micromegas(G4LogicalVolume*);

  //! add geometry node
  /*! this handles the internal (module/strips) segmentation, needed for tracking*/
  void add_geometry_node();

  //! detector parameters
  PHParameters* m_Params = nullptr;

  //! active volumes, and mapping to layer
  /*! it is needed in the stepping action to map a volume to a given layer */
  std::map<G4VPhysicalVolume*, int> m_activeVolumes;

  //! also store passive volumes
  std::set<G4VPhysicalVolume*> m_passiveVolumes;
  
  //! super detector name
  std::string m_SuperDetector;

  //! first layer number
  /* there are two layers in the detector */
  int m_FirstLayer = 0;

};

#endif // MICROMEGASDETECTOR_H
