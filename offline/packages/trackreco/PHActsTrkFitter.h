/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawleyu@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"

#include <trackbase/TrkrDefs.h>

#include <Acts/Geometry/GeometryContext.hpp>

#include <TMatrixDfwd.h>                      // for TMatrixD

#include <map>
#include <memory>                // for shared_ptr
#include <string>
#include <vector>


class SvtxTrackMap;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;
class TrkrClusterContainer;
class TGeoManager;
class TGeoNode;

namespace FW {
  class IBaseDetector;
  class IContextDecorator;
}

namespace Acts {
  class Surface;
}

//! \brief		Refit SvtxTracks with Acts.
class PHActsTrkFitter : public PHTrackFitting
{
 public:
  /*!
	 */

  //! Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  //! dtor
  ~PHActsTrkFitter();

  //!Initialization, called for initialization
  //int Init(PHCompositeNode*);

  //!Initialization Run, called for initialization of a run
  //  int InitRun(PHCompositeNode*);

  //!Process Event, called for each event
  //int process_event(PHCompositeNode*);

  //!End, write and close files
  int End(PHCompositeNode*);

int Setup(PHCompositeNode* topNode);

int Process();

  //Flags of different kinds of outputs
  enum Flag
  {
    //all disabled
    NONE = 0,
  };


 private:
  //! Event counter
  int _event;

  //! Get all the nodes
  int GetNodes(PHCompositeNode*);

  //!Create New nodes
  int CreateNodes(PHCompositeNode*);
  void BuildSiliconLayers();
void BuildTpcSurfaceMap();
  void isActive(TGeoNode *gnode);
  void MakeTGeoNodeMap(PHCompositeNode*);
  void getInttKeyFromNode(TGeoNode *gnode);
  void getMvtxKeyFromNode(TGeoNode *gnode);
  int MakeActsGeometry(int argc, char* argv[], FW::IBaseDetector& detector);
  TMatrixD GetMvtxCovarLocal(const unsigned int layer, const unsigned int staveid, const unsigned int chipid, TMatrixD world_err);
  TMatrixD GetInttCovarLocal(const unsigned int layer, const unsigned int staveid, const unsigned int chipid, TMatrixD world_err);
TMatrixD TransformCovarToLocal(const double ladderphi, TMatrixD world_err);
  TrkrDefs::hitsetkey GetMvtxHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world);
  TrkrDefs::hitsetkey GetInttHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world);

  /*
	 * fit track with SvtxTrack as input seed.
	 * \param intrack Input SvtxTrack
	 * \param invertex Input Vertex, if fit track as a primary vertex
	 */
  PHG4CylinderGeomContainer* _geom_container_mvtx;
  PHG4CylinderGeomContainer* _geom_container_intt;
  PHG4CylinderCellGeomContainer* _geom_container_tpc;

  SvtxTrackMap* _trackmap;
  TrkrClusterContainer* _clustermap;

  TGeoManager* _geomanager;

  Acts::GeometryContext  geo_ctxt;

  std::vector<std::shared_ptr<FW::IContextDecorator> > contextDecorators;

  std::map<TrkrDefs::hitsetkey, TGeoNode*> _cluster_node_map;
  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>> _cluster_surface_map;
  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>> _cluster_surface_map_tpc;

double SurfStepPhi;
double SurfStepZ;
double ModuleStepPhi;
double ModulePhiStart;

};

#endif
