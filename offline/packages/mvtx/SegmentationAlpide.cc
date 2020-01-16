/**
 * @file mvtx/SegmentationAlpide.cc
 * @author YCM, from Alice geom of 08/04/2019
 * @brief mvtx object with ALPIDE chip description
 */
#include "SegmentationAlpide.h"
#include <cstdio>

void SegmentationAlpide::print()
{
  printf("Pixel size: %.2f (along %d rows) %.2f (along %d columns) microns\n", PitchRow * 1e4, NRows, PitchCol * 1e4, NCols);
  printf("Passive edges: bottom: %.2f, top: %.2f, left/right: %.2f microns\n",
         PassiveEdgeReadOut * 1e4, PassiveEdgeTop * 1e4, PassiveEdgeSide * 1e4);
  printf("Active/Total size: %.6f/%.6f (rows) %.6f/%.6f (cols) cm\n", ActiveMatrixSizeRows, SensorSizeRows,
         ActiveMatrixSizeCols, SensorSizeCols);
}
