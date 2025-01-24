#ifndef G4MAIN_PHG4PROCESSMAPPHYSICS_H
#define G4MAIN_PHG4PROCESSMAPPHYSICS_H

#include "PHG4MCProcessDefs.h"

#include <map>
#include <string>

class G4VProcess;

class PHG4ProcessMapPhysics
{
 public:
  // explicit constructor
  PHG4ProcessMapPhysics();
  /// default copy constructor
  PHG4ProcessMapPhysics(const PHG4ProcessMapPhysics& rhs) = delete;
  /// default copy asigment
  PHG4ProcessMapPhysics& operator=(const PHG4ProcessMapPhysics& rhs) = delete;
  /// default move constructor
  PHG4ProcessMapPhysics(PHG4ProcessMapPhysics&& rhs) = delete;
  /// default move assigment constructor
  PHG4ProcessMapPhysics& operator=(PHG4ProcessMapPhysics&& rhs) = delete;
  /// default destructor
  ~PHG4ProcessMapPhysics() = default;

  // static access method
  static PHG4ProcessMapPhysics& Instance()
  {
    static PHG4ProcessMapPhysics fgInstance;
    return fgInstance;
  }

  // get methods
  PHG4MCProcess GetMCProcess(const G4VProcess* process) const;
  std::string_view GetMCProcessName(const G4VProcess* process) const;

 private:
  void FillMap();
};

#endif  // G4MAIN_PHG4PROCESSMAPPHYSICS_H
