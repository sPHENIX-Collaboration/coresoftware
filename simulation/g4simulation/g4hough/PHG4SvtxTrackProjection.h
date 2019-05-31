#ifndef G4HOUGH_PHG4SVTXTRACKPROJECTION_H
#define G4HOUGH_PHG4SVTXTRACKPROJECTION_H

//===========================================================
/// \file PHG4SvtxTrackProjection.h
/// \brief Projects into calorimeters and fills track cal fields
/// \author Mike McCumber
//===========================================================

#include "PHG4HoughTransform.h"

#include <trackbase_historic/SvtxTrack.h>

// PHENIX includes
#include <fun4all/SubsysReco.h>

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
  
  float get_mag_field() const          {return _magfield;}
  void  set_mag_field(float magfield) {_magfield = magfield;}
  
 private:

  PHG4HoughTransform _hough;
  int _num_cal_layers;
  std::vector<SvtxTrack::CAL_LAYER> _cal_types;
  std::vector<std::string> _cal_names;
  std::vector<float> _cal_radii;
  double _magfield;
};

#endif // G4HOUGH_PHG4INITZVERTEXING_H
