/*!
 *  \file		PHGenFitTrackProjection.h
 *  \brief		Projects into calorimeters and fills track cal fields using GenFit
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef G4HOUGH_PHGENFITTRACKPROJECTION_H
#define G4HOUGH_PHGENFITTRACKPROJECTION_H

#include <trackbase_historic/SvtxTrack.h>

#include <fun4all/SubsysReco.h>

// std includes
#include <vector>
#include <string>

// forward declarations
class PHCompositeNode;

namespace PHGenFit {
	class Fitter;
}

/// \class PHGenFitTrackProjection
///
/// \brief Projects into calorimeters and fills track cal fields
///
class PHGenFitTrackProjection : public SubsysReco
{

 public:
 
  PHGenFitTrackProjection(const std::string &name = "PHGenFitTrackProjection", const int pid_guess = 211);
  ~PHGenFitTrackProjection() override {}
		
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  int get_pid_guess() const {
    return _pid_guess;
  }
  
  void set_pid_guess(int pidGuess) {
    _pid_guess = pidGuess;
  }
  
  void use_poscalib_cemc_clusters(bool calib) {
    _use_poscalib_cemc = calib;
  }

 private:

  PHGenFit::Fitter * _fitter;
  int _pid_guess;

  bool _use_poscalib_cemc = false;

  int _num_cal_layers;
  std::vector<SvtxTrack::CAL_LAYER> _cal_types;
  std::vector<std::string> _cal_names;
  std::vector<float> _cal_radii;
};

#endif // __PHGENFITTRACKPROJECTION_H__
