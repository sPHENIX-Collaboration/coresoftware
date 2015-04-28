#ifndef __ZHOUGH_GRFITTER__
#define __ZHOUGH_GRFITTER__


#include "ZHough.h"


class ZHough_Cylindrical : public ZHough
{
  public:
    ZHough_Cylindrical(unsigned int z0_nbin, unsigned int theta_nbin, ZResolution& min_resolution, ZResolution& max_resolution, ZRange& range, double sxy, double sz);
    ~ZHough_Cylindrical();
    
    void findTracks(unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, std::vector<SimpleTrack3D>& tracks, std::vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, std::vector<float>& params);
    
    void findTracksCombo_noVertex(std::vector<SimpleHit3D>& comb_hits,  unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, std::vector<SimpleTrack3D>& tracks, std::vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, std::vector<float>& params);
    
    void findTracksCombo_withVertex(std::vector<SimpleHit3D>& comb_hits,  unsigned int zoomlevel, float xydiffcut, unsigned int max_hits, unsigned int tracks_per_hit, double chi2_cut, float max_kappa_cut, std::vector<SimpleTrack3D>& tracks, std::vector<SimpleHit3D>& hits, float min_k, float max_k, float min_phi, float max_phi, float min_d, float max_d, float min_z0, float max_z0, float min_theta, float max_theta, std::vector<float>& params);
    
    void setNLayers(unsigned int n){nlayers=n;}
    
    double fitTrack(SimpleTrack3D& track, std::vector<double>& chi2_hit);
    
    void setLayerResolution(std::vector<double>& lxy, std::vector<double>& lz)
    {
      layer_xy_resolution = lxy;
      layer_z_resolution = lz;
    }
    
    void setVertexResolution(double vxy, double vz)
    {
      vertex_sigma_xy = vxy;
      vertex_sigma_z = vz;
    }
    
    double getLayerResolution_xy(unsigned int l){return layer_xy_resolution[l];}
    double getLayerResolution_z(unsigned int l){return layer_z_resolution[l];}
    
  private:
    unsigned int nlayers;
    std::vector<double> layer_xy_resolution;
    std::vector<double> layer_z_resolution;
    double vertex_sigma_xy, vertex_sigma_z;
};




#endif
