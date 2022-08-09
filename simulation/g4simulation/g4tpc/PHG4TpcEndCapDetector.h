// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4TPCENDCAPDETECTOR_H
#define PHG4TPCENDCAPDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>

#include <set>
#include <string>  // for string
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;
class G4AssemblyVolume;
class PHG4TpcEndCapDisplayAction;

class PHG4TpcEndCapDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4TpcEndCapDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4TpcEndCapDetector() override;

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
  PHParameters *m_Params = nullptr;
  PHG4TpcEndCapDisplayAction *m_DisplayAction = nullptr;

  // active volumes
  std::set<G4LogicalVolume *> m_LogicalVolumesSet;

  std::string m_SuperDetector;

  G4AssemblyVolume *m_EndCapAssembly = nullptr;

  G4AssemblyVolume *ConstructEndCapAssembly();

  void ConstructWagonWheel(G4AssemblyVolume *assmeblyvol,
                           G4double &z_start);  // careful z_start is modified and being used later

  void ConstructElectronics(G4AssemblyVolume *assmeblyvol,
                            G4double z_start);

  void
  AddLayer(  //
      G4AssemblyVolume *assmeblyvol,
      G4double &z_start,
      const std::string &_name,        //! name base for this layer
      std::string _material,           //! material name in G4
      G4double _depth,                 //! depth in G4 units
      double _percentage_filled = 100  //! percentage filled//
  );

  void CreateCompositeMaterial(               //
      std::string compositeName,              //! desired name for the new material
      std::vector<std::string> materialName,  //! vector of the names of the component materials in G4
      std::vector<double> thickness           //! thickness of this particular layer (assuming 100 percent filled)
  );
};

#endif  // PHG4TPCENDCAPDETECTOR_H
