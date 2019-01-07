/*!
 *  \file		PHG4GenFitTrackProjection.h
 *  \brief		Projects into calorimeters and fills track cal fields using GenFit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4GENFITTRACKPROJECTION_H__
#define __PHG4GENFITTRACKPROJECTION_H__

#include "PHG4HoughTransform.h"

#include <trackbase_historic/SvtxTrack.h>

// PHENIX includes
#include <fun4all/SubsysReco.h>

// std includes
#include <vector>
#include <string>

// forward declarations
class PHCompositeNode;

namespace PHGenFit {
	class Fitter;
}

/// \class PHG4GenFitTrackProjection
///
/// \brief Projects into calorimeters and fills track cal fields
///
class PHG4GenFitTrackProjection : public SubsysReco
{

 public:
 
  PHG4GenFitTrackProjection(const std::string &name = "PHG4GenFitTrackProjection", const int pid_guess = 211);
  virtual ~PHG4GenFitTrackProjection() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

	int get_pid_guess() const {
		return _pid_guess;
	}

	void set_pid_guess(int pidGuess) {
		_pid_guess = pidGuess;
	}


 private:

  PHGenFit::Fitter * _fitter;
  int _pid_guess;

  int _num_cal_layers;
  std::vector<SvtxTrack::CAL_LAYER> _cal_types;
  std::vector<std::string> _cal_names;
  std::vector<float> _cal_radii;
};

#endif // __PHG4GENFITTRACKPROJECTION_H__
