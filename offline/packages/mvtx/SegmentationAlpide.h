/**
 * @file mvtx/SegmentationAlpide.h
 * @author YCM, from ALICE geom of 08/04/2019
 * @brief mvtx object with ALPIDE chip description
 */
#ifndef MVTX_ALPIDE_SEGMENTATION_H
#define MVTX_ALPIDE_SEGMENTATION_H

#include <TVector3.h>
#include <iostream>

class SegmentationAlpide
{
 public:
  static constexpr int   NCols = 1024;
  static constexpr int   NRows =  512;
  static constexpr int   NPixels = NRows * NCols;
  static constexpr float PitchCol = 29.24e-4;
  static constexpr float PitchRow = 26.88e-4;
  static constexpr float PassiveEdgeReadOut = 0.12f;              // width of the readout edge (Passive bottom)
  static constexpr float PassiveEdgeTop = 37.44e-4;               // Passive area on top
  static constexpr float PassiveEdgeSide = 29.12e-4;              // width of Passive area on left/right of the sensor
  static constexpr float ActiveMatrixSizeCols = PitchCol * NCols; // Active size along columns
  static constexpr float ActiveMatrixSizeRows = PitchRow * NRows; // Active size along rows

  // effective thickness of sensitive layer, accounting for charge collection non-unifoemity, https://alice.its.cern.ch/jira/browse/AOC-46
  static constexpr float SensorLayerThicknessEff = 22.e-4;
  static constexpr float SensorLayerThickness = 30.e-4;                                               // effective thickness of sensitive part
  static constexpr float SensorSizeCols = ActiveMatrixSizeCols + PassiveEdgeSide + PassiveEdgeSide;   // SensorSize along columns
  static constexpr float SensorSizeRows = ActiveMatrixSizeRows + PassiveEdgeTop + PassiveEdgeReadOut; // SensorSize along rows

  SegmentationAlpide() = default;
  ~SegmentationAlpide() = default;

  /// Transformation from Geant detector centered local coordinates (cm) to
  /// Pixel cell numbers iRow and iCol.
  /// Returns kTRUE if point x,z is inside sensitive volume, kFALSE otherwise.
  /// A value of -1 for iRow or iCol indicates that this point is outside of the
  /// detector segmentation as defined.
  /// @param float x Detector local coordinate x in cm with respect to
  /// the center of the sensitive volume.
  /// @param float z Detector local coordinate z in cm with respect to
  /// the center of the sensitive volulme.
  /// @param int iRow Detector x cell coordinate. Has the range 0 <= iRow < mNumberOfRows
  /// @param int iCol Detector z cell coordinate. Has the range 0 <= iCol < mNumberOfColumns
  static bool localToDetector(float x, float z, int& iRow, int& iCol);
  /// same but w/o check for row/column range
  static void localToDetectorUnchecked(float xRow, float zCol, int& iRow, int& iCol);

  /// Transformation from Detector cell coordiantes to Geant detector centered
  /// local coordinates (cm)
  /// @param int iRow Detector x cell coordinate. Has the range 0 <= iRow < mNumberOfRows
  /// @param int iCol Detector z cell coordinate. Has the range 0 <= iCol < mNumberOfColumns
  /// @param float x Detector local coordinate x in cm with respect to the
  /// center of the sensitive volume.
  /// @param float z Detector local coordinate z in cm with respect to the
  /// center of the sensitive volulme.
  /// If iRow and or iCol is outside of the segmentation range a value of -0.5*Dx()
  /// or -0.5*Dz() is returned.
  static bool detectorToLocal(int iRow, int iCol, float& xRow, float& zCol);
  static bool detectorToLocal(float row, float col, float& xRow, float& zCol);
  static bool detectorToLocal(float row, float col, TVector3& loc);

  // same but w/o check for row/col range
  static void detectorToLocalUnchecked(int iRow, int iCol, float& xRow, float& zCol);
  static void detectorToLocalUnchecked(float row, float col, float& xRow, float& zCol);
  static void detectorToLocalUnchecked(float row, float col, TVector3& loc);

  static constexpr float getFirstRowCoordinate()
  {
    return 0.5 * ((ActiveMatrixSizeRows - PassiveEdgeTop + PassiveEdgeReadOut) - PitchRow);
  }
  static constexpr float getFirstColCoordinate() { return 0.5 * (PitchCol - ActiveMatrixSizeCols); }

  static void print();
};

//_________________________________________________________________________________________________
inline void SegmentationAlpide::localToDetectorUnchecked(float xRow, float zCol, int& iRow, int& iCol)
{
  // convert to row/col w/o over/underflow check
  xRow = 0.5 * (ActiveMatrixSizeRows - PassiveEdgeTop + PassiveEdgeReadOut) - xRow; // coordinate wrt top edge of Active matrix
  zCol += 0.5 * ActiveMatrixSizeCols;                                               // coordinate wrt left edge of Active matrix
  iRow = int(xRow / PitchRow);
  iCol = int(zCol / PitchCol);
  if (xRow < 0)
    iRow -= 1;
  if (zCol < 0)
    iCol -= 1;
}

//_________________________________________________________________________________________________
inline bool SegmentationAlpide::localToDetector(float xRow, float zCol, int& iRow, int& iCol)
{
  // convert to row/col
  xRow = 0.5 * (ActiveMatrixSizeRows - PassiveEdgeTop + PassiveEdgeReadOut) - xRow; // coordinate wrt left edge of Active matrix
  zCol += 0.5 * ActiveMatrixSizeCols;                                               // coordinate wrt bottom edge of Active matrix
  if (xRow < 0 || xRow >= ActiveMatrixSizeRows || zCol < 0 || zCol >= ActiveMatrixSizeCols) {
    iRow = iCol = -1;
    return false;
  }
  iRow = int(xRow / PitchRow);
  iCol = int(zCol / PitchCol);
  return true;
}

//_________________________________________________________________________________________________
inline void SegmentationAlpide::detectorToLocalUnchecked(int iRow, int iCol, float& xRow, float& zCol)
{
  xRow = getFirstRowCoordinate() - iRow * PitchRow;
  zCol = iCol*PitchCol + getFirstColCoordinate();
}

//_________________________________________________________________________________________________
inline void SegmentationAlpide::detectorToLocalUnchecked(float row, float col, float& xRow, float& zCol)
{
  xRow = getFirstRowCoordinate() - row * PitchRow;
  zCol = col * PitchCol + getFirstColCoordinate();
}

//_________________________________________________________________________________________________
inline void SegmentationAlpide::detectorToLocalUnchecked(float row, float col, TVector3& loc)
{
  loc.SetXYZ(getFirstRowCoordinate() - row * PitchRow, 0.f, col * PitchCol + getFirstColCoordinate());
}

//_________________________________________________________________________________________________
inline bool SegmentationAlpide::detectorToLocal(int iRow, int iCol, float& xRow, float& zCol)
{
  if (iRow < 0 || iRow >= NRows || iCol<0 || iCol >= NCols)
    return false;
  detectorToLocalUnchecked(iRow,iCol,xRow,zCol);
  return true;
}

//_________________________________________________________________________________________________
inline bool SegmentationAlpide::detectorToLocal(float row, float col, float& xRow, float& zCol)
{
  if (row < 0 || row >= NRows || col < 0 || col >= NCols)
    return false;
  detectorToLocalUnchecked(row, col, xRow, zCol);
  return true;
}

//_________________________________________________________________________________________________
inline bool SegmentationAlpide::detectorToLocal(float row, float col, TVector3& loc)
{
  if (row < 0 || row >= NRows || col < 0 || col >= NCols)
    return false;
  detectorToLocalUnchecked(row, col, loc);
  return true;
}

#endif
