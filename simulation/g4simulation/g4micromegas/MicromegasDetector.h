// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MICROMEGASDETECTOR_H
#define MICROMEGASDETECTOR_H

/*!
 * \file MicromegasDetector.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4main/PHG4Detector.h>

#include <set>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

class MicromegasDetector : public PHG4Detector
{
 public:

  //! constructor
  MicromegasDetector(PHG4Subsystem*, PHCompositeNode*, PHParameters*, const std::string&);

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  //! set layer number
  void set_layer( int layer ) { m_layer = layer; }
  
  //! get layer number
  int get_layer() const { return m_layer; }
  
  private:
  
  // create needed material
  void create_materials() const;
  
  PHParameters *m_Params = nullptr;

  // active volumes
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;

  //! detector layer number
  int m_layer = 0;
  
};

#endif // MICROMEGASDETECTOR_H
