/**
 * @file mvtx/SegmentationAlpide.cc
 * @author YCM, from Alice geom of 08/04/2019
 * @brief mvtx object with ALPIDE chip description
 */
#include "SegmentationAlpide.h"

#include <boost/format.hpp>

#include <iostream>

void SegmentationAlpide::print()
{
  std::cout << (boost::format("Pixel size: %.2f (along %d rows) %.2f (along %d columns) microns") % (PitchRow * 1e4) % NRows % (PitchCol * 1e4) % NCols).str()
            << std::endl;
  std::cout << (boost::format("Passive edges: bottom: %.2f, top: %.2f, left/right: %.2f microns") % (PassiveEdgeReadOut * 1e4) % (PassiveEdgeTop * 1e4) % (PassiveEdgeSide * 1e4)).str()
            << std::endl;
  std::cout << (boost::format("Active/Total size: %.6f/%.6f (rows) %.6f/%.6f (cols) cm") % ActiveMatrixSizeRows % SensorSizeRows % ActiveMatrixSizeCols % SensorSizeCols).str()
            << std::endl;
}
