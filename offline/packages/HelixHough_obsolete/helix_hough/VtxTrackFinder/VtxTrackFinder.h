#ifndef __VTXTRACKFINDER__
#define __VTXTRACKFINDER__

#include "HelixHough.h"
#include <vector>
#include <set>

class VtxTrackFinder : public HelixHough
{
 public:
  VtxTrackFinder(std::vector<float>& detrad, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range);
  virtual ~VtxTrackFinder(){}
  
  void finalize(std::vector<SimpleTrack3D>& input, std::vector<SimpleTrack3D>& output);
  void findTracks(std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks);
  void initEvent(std::vector<SimpleHit3D>& hits)
  {
    combos.clear();
  }
  
  void setUsingVertex(bool usevtx){using_vertex = usevtx;}
  
  void setLayerResolution(std::vector<float>& lxy, std::vector<float>& lz)
  {
    layer_xy_resolution = lxy;
    layer_z_resolution = lz;
  }
    
  void setVertexResolution(float vxy, float vz)
  {
    vertex_sigma_xy = vxy;
    vertex_sigma_z = vz;
  }
  
  double fitTrack(SimpleTrack3D& track, std::vector<double>& chi2_hit);
  
  void setChi2Cut(double c){chi2_cut=c;}
  
 private:
  bool using_vertex;
  std::vector<float> detector_radii;
  std::vector<float> layer_xy_resolution;
  std::vector<float> layer_z_resolution;
  double vertex_sigma_xy, vertex_sigma_z;
  double chi2_cut;
  std::set<std::vector<unsigned int> > combos; 
};

#endif
