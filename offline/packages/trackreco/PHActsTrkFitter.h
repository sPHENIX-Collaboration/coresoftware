/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawleyu@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"
#include <fun4all/SubsysReco.h>

#include <cstddef>              // for NULL
#include <memory>                // for shared_ptr
#include <string>
#include <vector>
#include <map>

class TClonesArray;

class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4CylinderGeomContainer;
class TrkrClusterContainer;
class TTree;
class TGeoManager;
class TGeoNode;

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

  void isActive(TGeoNode *gnode);

  /*
	 * fit track with SvtxTrack as input seed.
	 * \param intrack Input SvtxTrack
	 * \param invertex Input Vertex, if fit track as a primary vertex
	 */
  PHG4CylinderGeomContainer* _geom_container_mvtx;
  PHG4CylinderGeomContainer* _geom_container_intt;
  SvtxTrackMap* _trackmap;
  TrkrClusterContainer* _clustermap;

  TGeoManager* _geomanager;
};

#endif
