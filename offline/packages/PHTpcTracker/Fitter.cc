/*!
 *  \file       Fitter.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

//PHGenFit2
#include "Fitter.h"


//GenFit
#include <GenFit/Exception.h>
#include <GenFit/FieldManager.h>
#include <GenFit/KalmanFitter.h>
#include <GenFit/MaterialEffects.h>
#include <GenFit/TGeoMaterialInterface.h>

//GenFitExp
#include <genfitexp/Field.h>

namespace PHGenFit2
{
  Fitter::Fitter(TGeoManager* tgeo_manager, PHField* field)
  //    : _tgeo_manager(tgeo_manager)
  //    , _field(field)
  {
    genfit::Exception::quiet(true);
    genfit::Field* fieldMap = new genfit::Field(field);
    genfit::FieldManager::getInstance()->init(fieldMap);

    // TODO: vacuum option for material effects?
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

    // init fitters
    _fitter5 = new genfit::KalmanFitter(5 /* maxIterations */);
    _fitter1 = new genfit::KalmanFitter(1 /* maxIterations */);
  }

  Fitter::~Fitter()
  {
    delete _fitter5;
    delete _fitter1;
  }

}  // namespace PHGenFit2
