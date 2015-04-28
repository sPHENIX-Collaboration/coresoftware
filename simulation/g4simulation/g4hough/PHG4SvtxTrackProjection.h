#ifndef __PHG4SVTXTRACKPROJECTION_H__
#define __PHG4SVTXTRACKPROJECTION_H__

//===========================================================
/// \file PHG4SvtxTrackProjection.h
/// \brief Projects into calorimeters and fills track cal fields
/// \author Mike McCumber
//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>

// PHG4 includes
#include <PHG4HoughTransform.h>

// std includes
#include <vector>
#include <string>

// forward declarations
class PHCompositeNode;

/// \class PHG4SvtxTrackProjection
///
/// \brief Projects into calorimeters and fills track cal fields
///
class PHG4SvtxTrackProjection : public SubsysReco
{

 public:
 
  PHG4SvtxTrackProjection(const std::string &name = "PHG4SvtxTrackProjection");
  virtual ~PHG4SvtxTrackProjection() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:

  PHG4HoughTransform _hough;
  int _num_cal_layers;
  std::vector<std::string> _cal_names;
  std::vector<float> _cal_radii;
  double _mag_extent;
};

#endif // __PHG4SVTXTRACKPROJECTION_H__
