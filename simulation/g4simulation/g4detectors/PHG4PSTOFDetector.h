#ifndef PHG4PSTOFDetector_h
#define PHG4PSTOFDetector_h

#include <g4main/PHG4Detector.h>

// cannot fwd declare G4RotationMatrix, it is a typedef pointing to clhep
#include <Geant4/G4RotationMatrix.hh>

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>

#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHG4Parameters;

class PHG4PSTOFDetector: public PHG4Detector
{
//  typedef CGAL::Exact_circular_kernel_2             Circular_k;
//  typedef CGAL::Point_2<Circular_k>                 Point_2;

  public:

  //! constructor
  PHG4PSTOFDetector( PHCompositeNode *Node, PHG4Parameters *params, const std::string &dnam="PSTOF");

  //! destructor
  virtual ~PHG4PSTOFDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  //int IsInPSTOF(G4VPhysicalVolume*) const;
  int IsInPSTOF(G4LogicalVolume*) const;
  //@}

  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  //int get_Layer() const {return layer;}

  protected:

  int active;
  G4LogicalVolume *active_volume;

  //int layer;
  //std::vector<G4VSolid *> scinti_tiles_vec; 
  //std::set<G4VPhysicalVolume *>steel_absorber_vec;

  std::string superdetector;
};

#endif
