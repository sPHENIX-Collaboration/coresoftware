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
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class PHG4MicromegasDetector : public PHG4Detector
{
 public:

  //* constructor
  PHG4MicromegasDetector(PHG4Subsystem*, PHCompositeNode*, PHParameters*, const std::string&);

  //* construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //* set first layer number
  void set_first_layer( int layer ) { m_first_layer = layer; }

  //* get first layer number
  int get_first_layer() const { return m_first_layer; }

  //* returns true if passed volume is an active volume of this detector
  bool IsInDetector(G4VPhysicalVolume*) const;

  //* return layer associated to a given volume, or -1 if invalid
  int get_layer(G4VPhysicalVolume*) const;

  //* super detector name
  void SuperDetector(const std::string &name) { m_superdetector = name; }

  //* super detector name
  const std::string SuperDetector() const { return m_superdetector; }

  private:

  //* create needed material
  void create_materials() const;

  //* detector parameters
  PHParameters* m_Params = nullptr;

  //* active volumes, and mapping to layer
  /** it is needed in the stepping action to map a volume to a given layer */
  std::map<G4VPhysicalVolume*, int> m_PhysicalVolumes;

  //* super detector name
  std::string m_superdetector;

  //* first layer number
  /* there are two layers in the detector */
  int m_first_layer = 0;

};

#endif // MICROMEGASDETECTOR_H
