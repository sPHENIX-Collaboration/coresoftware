/*!
 *  \file		Fitter.h
 *  \brief		Fitter class handles setups for the fitting.
 *  \details	Fitter class handles setups for the fitting like Geometry, Fields, fitter choice, etc.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_FITTER_H
#define PHGENFIT_FITTER_H

// needed, it crashes on Ubuntu using singularity with local cvmfs install
// shared pointer later on uses this, forward declaration does not cut it
#include <phgenfit/Track.h> 

#include <GenFit/EventDisplay.h>
#include "GenFit/Exception.h"

#include <string>

//BOOST
//#include<boost/make_shared.hpp>
//
//#define SMART(expr) boost::shared_ptr<expr>
//#define NEW(expr) boost::make_shared<expr>

class TGeoManager;
class PHField;

namespace genfit
{
class AbsKalmanFitter;
class AbsBField;
}  // namespace genfit

namespace PHGenFit
{

class Fitter
{
 public:
  enum FitterType
  {
    KalmanFitter,
    KalmanFitterRefTrack,
    DafSimple,
    DafRef
  };
  enum TrackRepType
  {
    RKTrackRep
  };

  //! Default constructor
  Fitter(const std::string& tgeo_file_name,
         const PHField* field,
         const std::string& fitter_choice = "KalmanFitterRefTrack",
         const std::string& track_rep_choice = "RKTrackRep",
         const bool doEventDisplay = false);

  Fitter(TGeoManager* tgeo_manager,
         genfit::AbsBField* fieldMap,
         const std::string& fitter_choice = "KalmanFitterRefTrack",
         const std::string& track_rep_choice = "RKTrackRep",
         const bool doEventDisplay = false);

  Fitter(TGeoManager* tgeo_manager,
         genfit::AbsBField* fieldMap,
         const PHGenFit::Fitter::FitterType& fitter_choice = PHGenFit::Fitter::KalmanFitter,
         const PHGenFit::Fitter::TrackRepType& track_rep_choice = PHGenFit::Fitter::RKTrackRep,
         const bool doEventDisplay = false);

  //! Default destructor
  ~Fitter();

  static Fitter* getInstance(const std::string& tgeo_file_name,
                             const PHField* field,
                             const std::string& fitter_choice = "KalmanFitterRefTrack",
                             const std::string& track_rep_choice = "RKTrackRep",
                             const bool doEventDisplay = false);

  static Fitter* getInstance(TGeoManager* tgeo_manager,
                             const PHField* field,
                             const std::string& fitter_choice = "KalmanFitterRefTrack",
                             const std::string& track_rep_choice = "RKTrackRep",
                             const bool doEventDisplay = false);

  static Fitter* getInstance(TGeoManager* tgeo_manager,
                             const PHField* field,
                             const PHGenFit::Fitter::FitterType& fitter_choice = PHGenFit::Fitter::KalmanFitter,
                             const PHGenFit::Fitter::TrackRepType& track_rep_choice = PHGenFit::Fitter::RKTrackRep,
                             const bool doEventDisplay = false);

  int processTrack(PHGenFit::Track* track, const bool save_to_evt_disp = false);

  int displayEvent();

  bool is_do_Event_Display() const
  {
    return _doEventDisplay;
  }

  void set_do_Event_Display(bool doEventDisplay)
  {
    _doEventDisplay = doEventDisplay;
    if (!_display && _doEventDisplay)
      _display = genfit::EventDisplay::getInstance();
  }

  genfit::EventDisplay* getEventDisplay()
  {
    return _display;
  }

  int get_verbosity() const
  {
    return verbosity;
  }

  void set_verbosity(int v)
  {
    this->verbosity = v;
    if (verbosity >= 1)
      genfit::Exception::quiet(false);
    else
      genfit::Exception::quiet(true);
  }

 private:
  /*!
	 * Verbose control:
	 * -1: Silient
	 * 0: Minimum
	 * 1: Errors only
	 * 2: Errors and Warnings
	 * 3: Verbose mode, long term debugging
	 */
  int verbosity;

  TGeoManager* _tgeo_manager;

  bool _doEventDisplay;

  genfit::EventDisplay* _display;
  genfit::AbsKalmanFitter* _fitter;

};  //class Fitter

}  // namespace PHGenFit

#endif
