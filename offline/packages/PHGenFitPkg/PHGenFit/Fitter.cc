/*!
 *  \file		Fitter.cc
 *  \brief		Fitter class handles setups for the fitting.
 *  \details	Fitter class handles setups for the fitting like Geometry, Fields, fitter choice, etc.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

//PHGenFit
#include "Fitter.h"

#include "Track.h"

//ROOT
#include <TGeoManager.h>
#include <RVersion.h>                      // for ROOT_VERSION, ROOT_VERSION...

//GenFit
#include <GenFit/AbsKalmanFitter.h>
#include <GenFit/DAF.h>
#include <GenFit/EventDisplay.h>
#include <GenFit/FieldManager.h>
#include <GenFit/FitStatus.h>
#include <GenFit/KalmanFitter.h>
#include <GenFit/KalmanFitterRefTrack.h>
#include <GenFit/MaterialEffects.h>
#include <GenFit/TGeoMaterialInterface.h>
#include <GenFit/Track.h>

//GenFitExp
#include <genfitexp/Field.h>

#include <cassert>
#include <cstddef>
#include <iostream>

namespace genfit { class AbsTrackRep; }

#define LogDEBUG(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogERROR(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWARNING(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

namespace PHGenFit
{
Fitter::Fitter(
    const std::string& tgeo_file_name,
    const PHField* field,
    const std::string& fitter_choice,
    const std::string& /*track_rep_choice*/,
    const bool doEventDisplay)
  : verbosity(1000)
  , _doEventDisplay(doEventDisplay)
{
  _tgeo_manager = new TGeoManager("Default", "Geane geometry");
  TGeoManager::Import(tgeo_file_name.data());

  assert(field);
  genfit::Field* fieldMap = new genfit::Field(field);

  genfit::FieldManager::getInstance()->init(
      fieldMap);
  genfit::MaterialEffects::getInstance()->init(
      new genfit::TGeoMaterialInterface());

  // init event display
  if (_doEventDisplay)
    _display = genfit::EventDisplay::getInstance();
  else
    _display = NULL;

  // init fitter
  if (fitter_choice.compare("KalmanFitterRefTrack") == 0)
    _fitter = new genfit::KalmanFitterRefTrack();
  else if (fitter_choice.compare("KalmanFitter") == 0)
    _fitter = new genfit::KalmanFitter();
  else if (fitter_choice.compare("DafSimple") == 0)
    _fitter = new genfit::DAF(false);
  else if (fitter_choice.compare("DafRef") == 0)
    _fitter = new genfit::DAF(true);
  else
    _fitter = new genfit::KalmanFitter();

  genfit::Exception::quiet(true);
}

Fitter::~Fitter()
{
  delete _fitter;
  //delete _tgeo_manager;
  //_tgeo_manager->Delete();
  delete _display;
}

int Fitter::processTrack(PHGenFit::Track* track, const bool save_to_evt_disp)
{
  genfit::Track* fitTrack = track->getGenFitTrack();

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)
  try
  {
    fitTrack->checkConsistency();
  }
  catch (genfit::Exception& e)
  {
    if (verbosity >= 2)
    {
      std::cerr << "genfit::Track::checkConsistency() failed!" << std::endl;
      std::cerr << e.what();
    }
    return -1;
  }
#else
  if (!fitTrack->checkConsistency())
  {
    if (verbosity >= 2) LogWARNING("genfit::Track::checkConsistency() failed!");
    return -1;
  }
#endif
  try
  {
    _fitter->processTrack(fitTrack);
  }
  catch (genfit::Exception& e)
  {
    if (verbosity >= 1)
    {
      std::cerr << "PHGenFit::Fitter::processTrack::Exception: \n";
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
    }
    return -1;
  }
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)
  try
  {
    fitTrack->checkConsistency();
  }
  catch (genfit::Exception& e)
  {
    if (verbosity >= 2)
    {
      std::cerr << "genfit::Track::checkConsistency() failed!" << std::endl;
      std::cerr << e.what();
    }
    return -1;
  }
#else

  if (!fitTrack->checkConsistency())
  {
    if (verbosity >= 2) LogWARNING("genfit::Track::checkConsistency() failed!");
    return -1;
  }
#endif
  genfit::AbsTrackRep* rep = fitTrack->getCardinalRep();
  if (!fitTrack->getFitStatus(rep)->isFitConverged())
  {
    if (verbosity >= 2) LogWARNING("Track could not be fitted successfully! Fit is not converged!");
    return -1;
  }

  if (_display and save_to_evt_disp)
    _display->addEvent(track->getGenFitTrack());

  return 0;
}

Fitter* Fitter::getInstance(const std::string& tgeo_file_name,
                            const PHField* field,
                            const std::string& fitter_choice, const std::string& track_rep_choice,
                            const bool doEventDisplay)
{
  TGeoManager* tgeo_manager = TGeoManager::Import(tgeo_file_name.data(), "Default");
  if (!tgeo_manager)
  {
    LogERROR("No TGeoManager found!");
    return NULL;
  }

  assert(field);
  genfit::Field* fieldMap = new genfit::Field(field);
  return new Fitter(tgeo_manager, fieldMap, fitter_choice, track_rep_choice, doEventDisplay);
}

Fitter::Fitter(TGeoManager* tgeo_manager, genfit::AbsBField* fieldMap,
               const PHGenFit::Fitter::FitterType& fitter_choice,
               const PHGenFit::Fitter::TrackRepType& /*track_rep_choice*/,
               const bool doEventDisplay)
  : verbosity(0)
  , _tgeo_manager(tgeo_manager)
  , _doEventDisplay(doEventDisplay)
{
  genfit::FieldManager::getInstance()->init(
      fieldMap);
  genfit::MaterialEffects::getInstance()->init(
      new genfit::TGeoMaterialInterface());

  // init event display
  if (_doEventDisplay)
    _display = genfit::EventDisplay::getInstance();
  else
    _display = NULL;

  // init fitter
  if (fitter_choice == PHGenFit::Fitter::KalmanFitter)
    _fitter = new genfit::KalmanFitter();
  else if (fitter_choice == PHGenFit::Fitter::KalmanFitterRefTrack)
    _fitter = new genfit::KalmanFitterRefTrack();
  if (fitter_choice == PHGenFit::Fitter::DafSimple)
    _fitter = new genfit::DAF(false);
  else if (fitter_choice == PHGenFit::Fitter::DafRef)
    _fitter = new genfit::DAF(true);
  else
  {
    _fitter = nullptr;
    LogERROR("This fitter not implemented!");
  }
}

Fitter* Fitter::getInstance(TGeoManager* tgeo_manager,
                            const PHField* field,
                            const std::string& fitter_choice, const std::string& track_rep_choice,
                            const bool doEventDisplay)
{
  if (!tgeo_manager)
  {
    LogERROR("No TGeoManager found!");
    return NULL;
  }

  assert(field);
  genfit::Field* fieldMap = new genfit::Field(field);
  return new Fitter(tgeo_manager, fieldMap, fitter_choice, track_rep_choice, doEventDisplay);
}

Fitter::Fitter(TGeoManager* tgeo_manager, genfit::AbsBField* fieldMap,
               const std::string& fitter_choice, const std::string& /*track_rep_choice*/,
               const bool doEventDisplay)
  : verbosity(0)
  , _tgeo_manager(tgeo_manager)
  , _doEventDisplay(doEventDisplay)
{
  genfit::FieldManager::getInstance()->init(
      fieldMap);
  genfit::MaterialEffects::getInstance()->init(
      new genfit::TGeoMaterialInterface());

  // init event display
  if (_doEventDisplay)
    _display = genfit::EventDisplay::getInstance();
  else
    _display = NULL;

  // init fitter
  if (fitter_choice.compare("KalmanFitterRefTrack") == 0)
    _fitter = new genfit::KalmanFitterRefTrack();
  else if (fitter_choice.compare("KalmanFitter") == 0)
    _fitter = new genfit::KalmanFitter();
  else if (fitter_choice.compare("DafSimple") == 0)
    _fitter = new genfit::DAF(false);
  else if (fitter_choice.compare("DafRef") == 0)
    _fitter = new genfit::DAF(true);
  else
  {
    _fitter = nullptr;
    LogERROR("This fitter not implemented!");
  }
}

int Fitter::displayEvent()
{
  if (_display)
    _display->open();
  else if (verbosity >= 0)
    LogERROR("No genfit::EventDisplay found!");

  return 0;
}

Fitter* Fitter::getInstance(TGeoManager* tgeo_manager,
                            const PHField* field,
                            const PHGenFit::Fitter::FitterType& fitter_choice,
                            const PHGenFit::Fitter::TrackRepType& track_rep_choice,
                            const bool doEventDisplay)
{
  if (!tgeo_manager)
  {
    LogERROR("No TGeoManager found!");
    return NULL;
  }

  assert(field);
  genfit::Field* fieldMap = new genfit::Field(field);

  return new Fitter(tgeo_manager, fieldMap, fitter_choice, track_rep_choice, doEventDisplay);
}

}  // namespace PHGenFit
