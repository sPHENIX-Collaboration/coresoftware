/*!
 *  \file		PHG4GenFitTrackProjection.h
 *  \brief		Projects into calorimeters and fills track cal fields using GenFit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHG4GENFITTRACKPROJECTION_H__
#define __PHG4GENFITTRACKPROJECTION_H__

#include "PHG4HoughTransform.h"

#include "SvtxTrack.h"

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
  
  float get_mag_field() const          {return _magfield;}
  void  set_mag_field(float magfield) {_magfield = magfield;}

	const std::string& get_mag_field_file_name() const {
		return _mag_field_file_name;
	}

	void set_mag_field_file_name(const std::string& magFieldFileName) {
		_mag_field_file_name = magFieldFileName;
	}

	float get_mag_field_re_scaling_factor() const {
		return _mag_field_re_scaling_factor;
	}

	void set_mag_field_re_scaling_factor(float magFieldReScalingFactor) {
		_mag_field_re_scaling_factor = magFieldReScalingFactor;
	}

	int get_pid_guess() const {
		return _pid_guess;
	}

	void set_pid_guess(int pidGuess) {
		_pid_guess = pidGuess;
	}

 private:

  PHGenFit::Fitter * _fitter;
  std::string _mag_field_file_name;
  float _mag_field_re_scaling_factor;

  int _pid_guess;

  PHG4HoughTransform _hough;
  int _num_cal_layers;
  std::vector<SvtxTrack::CAL_LAYER> _cal_types;
  std::vector<std::string> _cal_names;
  std::vector<float> _cal_radii;
  double _magfield;
  double _mag_extent;
};

#endif // __PHG4GENFITTRACKPROJECTION_H__
