#ifndef G4MAIN_PHG4PROCESSMAPPHYSICS_H
#define G4MAIN_PHG4PROCESSMAPPHYSICS_H

#include <cstdint>
#include <string>
#include <map>

class G4VProcess;

class PHG4ProcessMapPhysics {
 public:
  // explicit constructor
  PHG4ProcessMapPhysics();
  /// default copy constructor
  PHG4ProcessMapPhysics(const PHG4ProcessMapPhysics& rhs);
  /// default copy asigment
  PHG4ProcessMapPhysics& operator=(const PHG4ProcessMapPhysics& rhs) = delete;
  /// default move constructor
  PHG4ProcessMapPhysics(PHG4ProcessMapPhysics&& rhs) = delete;
  /// default move assigment constructor
  PHG4ProcessMapPhysics& operator=(PHG4ProcessMapPhysics&& rhs) = delete;
  /// default destructor
  ~PHG4ProcessMapPhysics() = default;

  PHG4ProcessMapPhysics& Instance();

 private:
  void FillMap();
};

#endif  // G4MAIN_PHG4PROCESSMAPPHYSICS_H
