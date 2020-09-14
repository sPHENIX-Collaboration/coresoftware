// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4TPCENDCAPDETECTOR_H
#define PHG4TPCENDCAPDETECTOR_H

#include <g4main/PHG4Detector.h>
#include <Geant4/G4Types.hh>

#include <set>
#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;
class G4AssemblyVolume;
class PHG4TpcEndCapDisplayAction;
class G4VSolid;

class PHG4TpcEndCapDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TpcEndCapDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4TpcEndCapDetector();

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

 private:
  PHParameters *m_Params;
  PHG4TpcEndCapDisplayAction *m_DisplayAction;

  // active volumes
  std::set<G4VPhysicalVolume *> m_PhysicalVolumesSet;

  std::string m_SuperDetector;

  G4AssemblyVolume *m_EndCapAssembly = nullptr;

  G4AssemblyVolume *ConstructEndCapAssembly();
  void ConstructWagonWheel(G4AssemblyVolume *assmeblyvol,
      G4double &z_start);

  void
  AddLayer(                            //
      G4AssemblyVolume * assmeblyvol,
      G4double & z_start,
      std::string _name,               //! name base for this layer
      std::string _material,           //! material name in G4
      G4double _depth,                   //! depth in G4 units
      double _percentage_filled = 100  //! percentage filled//
  );
};

#endif  // PHG4TPCENDCAPDETECTOR_H
