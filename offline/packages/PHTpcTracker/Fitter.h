/*!
 *  \file       Fitter.h
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#ifndef PHGENFIT2_FITTER_H
#define PHGENFIT2_FITTER_H

#include "Track.h"

#include <GenFit/KalmanFitter.h>
#include <GenFit/Track.h>

class TGeoManager;
class PHField;

namespace PHGenFit2
{
  class Fitter
  {
   public:
    Fitter(TGeoManager* tgeo_manager, PHField* field);
    ~Fitter();

    genfit::KalmanFitter* getFitter5() { return _fitter5; }
    genfit::KalmanFitter* getFitter1() { return _fitter5; }

    void processTrack5(PHGenFit2::Track* track) { return processTrack5(track->getGenFitTrack()); }
    void processTrack1(PHGenFit2::Track* track) { return processTrack1(track->getGenFitTrack()); }

    void processTrack5(genfit::Track* track) { return _fitter5->processTrack(track); }
    void processTrack1(genfit::Track* track) { return _fitter1->processTrack(track); }

    void processTrackPartially(PHGenFit2::Track* track, int startId = 0, int endId = -1)
    {
      return processTrackPartially(track->getGenFitTrack(), startId, endId);
    }
    void processTrackPartially(genfit::Track* track, int startId = 0, int endId = -1)
    {
      return _fitter1->processTrackPartially(track, track->getCardinalRep(), startId, endId);
    }

   private:
    //    TGeoManager* _tgeo_manager;
    //    PHField* _field;
    genfit::KalmanFitter* _fitter5;  // KalmanFitter with maxIterations = 5
    genfit::KalmanFitter* _fitter1;  // KalmanFitter with maxIterations = 1

  };  //class Fitter

}  // namespace PHGenFit2

#endif
