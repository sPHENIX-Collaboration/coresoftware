#include "PHG4MvtxMisalignment.h"

#include <array>
#include <cmath>   // for M_PI
#include <cstdio>  // for sprintf
#include <fstream>
#include <iostream>  // for operator<<, basic...
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>    // for tuple
#include <utility>  // for pair, make_pair
#include <vector>

PHG4MvtxMisalignment::PHG4MvtxMisalignment() = default;

std::vector<double> PHG4MvtxMisalignment::get_GlobalDisplacement()
{
  std::vector<double> globalDisplacement = {m_GlobalDisplacementX, m_GlobalDisplacementY, m_GlobalDisplacementZ};
  return globalDisplacement;
}

void PHG4MvtxMisalignment::LoadMvtxStaveAlignmentParameters()
{
  std::ifstream file(mvtxStaveAlignParamsFile);
  if (!file.is_open())
  {
    std::cout << "PHG4MvtxMisalignment::LoadMvtxStaveAlignmentParameters - ERROR - Could not open file " << mvtxStaveAlignParamsFile << std::endl;
    m_GlobalDisplacementX = 0.0;
    m_GlobalDisplacementY = 0.0;
    m_GlobalDisplacementZ = 0.0;
    return;
  }
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream iss(line);
    int layer;
    int stave;
    double alpha;
    double beta;
    double gamma;
    double x;
    double y;
    double z;
    if (!(iss >> layer >> stave >> alpha >> beta >> gamma >> x >> y >> z))
    {
      std::cout << "PHG4MvtxMisalignment::LoadMvtxStaveAlignmentParameters - ERROR - Could not read line " << line << std::endl;
      m_GlobalDisplacementX = 0.0;
      m_GlobalDisplacementY = 0.0;
      m_GlobalDisplacementZ = 0.0;
      continue;
    }

    m_GlobalDisplacementX += x;
    m_GlobalDisplacementY += y;
    m_GlobalDisplacementZ += z;
  }

  double Ntotstaves = 0.0;
  for (const auto &i : PHG4MvtxDefs::mvtxdat)
  {
    Ntotstaves += i[PHG4MvtxDefs::kNStave];
  }

  m_GlobalDisplacementX /= Ntotstaves;
  m_GlobalDisplacementY /= Ntotstaves;
  m_GlobalDisplacementZ /= Ntotstaves;

  file.close();
}
