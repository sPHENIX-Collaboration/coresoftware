/*!
 *  \file		PHRaveVertexing.h
 *  \brief		Refit SvtxTracks with PHGenFit.
 *  \details	Refit SvtxTracks with PHGenFit.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHRAVEVERTEXING_H
#define TRACKRECO_PHRAVEVERTEXING_H

#include <fun4all/SubsysReco.h>

#include <map>                   // for map, map<>::value_compare
#include <string>
#include <vector>

namespace genfit
{
class GFRaveVertex;
class GFRaveVertexFactory;
class Track;
} /* namespace genfit */

namespace PHGenFit
{
class Fitter;
} /* namespace PHGenFit */

class PHTimer;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class PHCompositeNode;
class PHG4TruthInfoContainer;

//! \brief		Refit SvtxTracks with PHGenFit.
class PHRaveVertexing : public SubsysReco
{
 public:
  typedef std::map<const genfit::Track*, unsigned int> GenFitTrackMap;

  //! Default constructor
  PHRaveVertexing(const std::string& name = "PHRaveVertexing");

  //! dtor
  ~PHRaveVertexing() override;

  //!Initialization, called for initialization
  int Init(PHCompositeNode*) override;

  //!Initialization Run, called for initialization of a run
  int InitRun(PHCompositeNode*) override;

  //!Process Event, called for each event
  int process_event(PHCompositeNode*) override;

  //!End, write and close files
  int End(PHCompositeNode*) override;

  const std::string& get_vertexing_method() const
  {
    return _vertexing_method;
  }

  void set_vertexing_method(const std::string& vertexingMethod)
  {
    _vertexing_method = vertexingMethod;
  }

  int get_primary_pid_guess() const
  {
    return _primary_pid_guess;
  }

  void set_primary_pid_guess(int primaryPidGuess)
  {
    _primary_pid_guess = primaryPidGuess;
  }

  bool is_over_write_svtxvertexmap() const
  {
    return _over_write_svtxvertexmap;
  }

  void set_over_write_svtxvertexmap(bool overWriteSvtxvertexmap)
  {
    _over_write_svtxvertexmap = overWriteSvtxvertexmap;
  }

  void set_svtxvertexmaprefit_node_name(const std::string & name) {_svtxvertexmaprefit_node_name = name;}

  double get_vertex_min_ndf() const
  {
    return _vertex_min_ndf;
  }

  void set_vertex_min_ndf(double vertexMinPT)
  {
    _vertex_min_ndf = vertexMinPT;
  }

  void set_nmvtx_clusters_required(unsigned int n)
  {
    _nmvtx_required = n;
  }

 private:
  //! Event counter
  int _event;

  //! Get all the nodes
  int GetNodes(PHCompositeNode*);

  //!Create New nodes
  int CreateNodes(PHCompositeNode*);

  genfit::Track* TranslateSvtxToGenFitTrack(SvtxTrack* svtx);

  //! Fill SvtxVertexMap from GFRaveVertexes and Tracks
  bool FillSvtxVertexMap(
      const std::vector<genfit::GFRaveVertex*>& rave_vertices,
      const GenFitTrackMap& gf_track_map);

  bool _over_write_svtxvertexmap;
  std::string _svtxvertexmaprefit_node_name;

  PHGenFit::Fitter* _fitter;

  int _primary_pid_guess;
  double _vertex_min_ndf;

  unsigned int _nmvtx_required = 3;  // require 3 or more mvtx clusters for track to be used in vertexing

  genfit::GFRaveVertexFactory* _vertex_finder;

  //! https://rave.hepforge.org/trac/wiki/RaveMethods
  std::string _vertexing_method;

  //! Input Node pointers
  SvtxTrackMap* _trackmap;
  SvtxVertexMap* _vertexmap;

  //! Output Node pointers
  SvtxVertexMap* _vertexmap_refit;

  PHTimer* _t_translate;
  PHTimer* _t_rave;
};

#endif
