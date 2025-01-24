#ifndef G4MAIN_PHG4PROCESSMAP_H
#define G4MAIN_PHG4PROCESSMAP_H

#include "PHG4MCProcessDefs.h"

#include <cstdint>
#include <map>
#include <string_view>

class G4VProcess;

class PHG4ProcessMap
{
 public:
  PHG4ProcessMap() = default;
  PHG4ProcessMap(const PHG4ProcessMap& rhs) = delete;
  PHG4ProcessMap(PHG4ProcessMap&& rhs) = delete;
  PHG4ProcessMap& operator=(const PHG4ProcessMap& rhs) = delete;
  PHG4ProcessMap& operator=(PHG4ProcessMap&& rhs) = delete;
  ~PHG4ProcessMap() { Clear(); }

  // static access method
  static PHG4ProcessMap& Instance()
  {
    static PHG4ProcessMap fgInstance;
    return fgInstance;
  };

  // methods
  bool Add(int subType, PHG4MCProcess mcProcess);
  void PrintAll() const;
  void Clear();

  // get methods
  PHG4MCProcess GetMCProcess(const G4VProcess* process) const;
  std::string_view GetMCProcessName(const G4VProcess* process) const;

 private:
  // methods
  bool IsDefined(int subType);

  // data members
  std::map<int, PHG4MCProcess> fMap;
};

#endif  // G4MAIN_PHG4PROCESSMAP_H
