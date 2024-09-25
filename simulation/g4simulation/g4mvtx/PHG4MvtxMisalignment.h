#ifndef G4MVTX_PHG4MVTXMISALIGNMENT_H
#define G4MVTX_PHG4MVTXMISALIGNMENT_H

#include "PHG4MvtxDefs.h"

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <utility>

class PHG4MvtxMisalignment
{
 public:
  PHG4MvtxMisalignment();
  virtual ~PHG4MvtxMisalignment() = default;

  std::vector<double> get_GlobalDisplacement();

 private:
  std::string mvtxStaveAlignParamsFile = "./MvtxStaveAlignmentParameters_Run2024.txt"; //TODO: either put this text file in CDB or create CDBTree)
  double m_GlobalDisplacementX = 0.;
  double m_GlobalDisplacementY = 0.;
  double m_GlobalDisplacementZ = 0.;

  void LoadMvtxStaveAlignmentParameters();

};

#endif