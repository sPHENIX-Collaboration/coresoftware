/*!
 *  \file       PHTpcVertexFinder.cc
 *  \brief      
 *  \author     Dmitry Arkhipkin <arkhipkin@gmail.com>
 */

#include "PHTpcVertexFinder.h"

#include "Track.h"  // for Track

#include <phool/PHLog.h>

#include <GenFit/GFRaveVertexFactory.h>  // for GFRaveVertexFactory

#include <TMatrixDSymfwd.h>  // for TMatrixDSym
#include <TMatrixTSym.h>     // for TMatrixTSym
#include <TVector3.h>        // for TVector3

#include <log4cpp/CategoryStream.hh>  // for CategoryStream

namespace genfit
{
  class GFRaveVertex;
}
namespace genfit
{
  class Track;
}

PHTpcVertexFinder::PHTpcVertexFinder()
  : _vertex_finder(nullptr)
{
}

PHTpcVertexFinder::~PHTpcVertexFinder()
{
  delete _vertex_finder;
}

std::vector<genfit::GFRaveVertex*> PHTpcVertexFinder::findVertices(std::vector<PHGenFit2::Track*>& gtracks)
{
  std::vector<genfit::GFRaveVertex*> rave_vertices;

  if (!_vertex_finder)
  {
    _vertex_finder = new genfit::GFRaveVertexFactory(1 /* verbosity? */);
    _vertex_finder->setMethod("avr-minweight:0.5-primcut:5-seccut:10");
    //		_vertex_finder->setMethod("avf-Tini:256-ratio:0.25-sigmacut:3");
    TMatrixDSym beam_spot_cov(3);
    beam_spot_cov(0, 0) = 0.2 * 0.2;      // dx * dx
    beam_spot_cov(1, 1) = 0.2 * 0.2;      // dy * dy
    beam_spot_cov(2, 2) = 100.0 * 100.0;  // dz * dz
    _vertex_finder->setBeamspot(TVector3(0, 0, 0), beam_spot_cov);
  }

  std::vector<genfit::Track*> gftracks;
  gftracks.reserve(gtracks.size());
  for (auto it = gtracks.begin(); it != gtracks.end(); ++it)
  {
    gftracks.push_back((*it)->getGenFitTrack());
  }

  LOG_DEBUG("tracking.PHTpcVertexFinder.findVertices") << "feeding " << gftracks.size() << " GenFit tracks into vertex finder";

  if (gftracks.size() < 2)
  {
    LOG_DEBUG("tracking.PHTpcVertexFinder.findVertices") << "less than two tracks, skipping vertexing";
    return rave_vertices;
  }

  try
  {
    _vertex_finder->findVertices(&rave_vertices, gftracks, true /* useBeamSpot */);
  }
  catch (...)
  {
    LOG_DEBUG("tracking.PHTpcVertexFinder.findVertices") << "Rave threw an exception, vertex finding failed";
    rave_vertices.clear();
  }

  /*
	for ( unsigned int i = 0, ilen = rave_vtx->getNTracks(); i++ ) {
		const genfit::Track* gftrack = rave_vtx->getParameters(i)->getTrack();
	}
	*/

  return rave_vertices;
}
