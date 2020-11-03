#ifndef G4TPC_DISTORTIONSTRUCT_H
#define G4TPC_DISTORTIONSTRUCT_H

#include <vector>

  class DistortionStruct
  {
    
    public:
    using List = std::vector<DistortionStruct>;
    
    /// constructor
    DistortionStruct() = default;
    
    float _r = 0;
    float _phi = 0;
    float _z = 0;
    
    float _dr = 0;
    float _dphi = 0;
    float _dz = 0;

  };

#endif  // G4TPC_DISTORTIONSTRUCT_H
