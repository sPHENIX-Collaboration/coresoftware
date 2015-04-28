#ifndef __FOURHITSEEDFINDER__
#define __FOURHITSEEDFINDER__

#include "HelixHough.h"
#include <vector>
#include <set>


class FourHitSeedFinder : public HelixHough
{
  public:
    FourHitSeedFinder(std::vector<float>& detrad, unsigned int n_phi, unsigned int n_d, unsigned int n_k, unsigned int n_dzdl, unsigned int n_z0, HelixResolution& min_resolution, HelixResolution& max_resolution, HelixRange& range);
    virtual ~FourHitSeedFinder(){}
    
    void finalize(std::vector<SimpleTrack3D>& input, std::vector<SimpleTrack3D>& output);
    void findTracks(std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks, const HelixRange& range);
    void initEvent(std::vector<SimpleHit3D>& hits, unsigned int min_hits)
    {
      combos.clear();
      combos_3_pass.clear();
      combos_3_fail.clear();
    }
    
    void findTracks_3_4(std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks, const HelixRange& range);
    void findTracks_6(std::vector<SimpleHit3D>& hits, std::vector<SimpleTrack3D>& tracks, const HelixRange& range);
    
    void setUsingVertex(bool usevtx){using_vertex = usevtx;}
    
    void setLayerResolution(std::vector<float>& lxy, std::vector<float>& lz)
    {
      layer_xy_resolution = lxy;
      layer_z_resolution = lz;
    }
    
    // for each detector layer, the material budget in radiation lengths
    void setLayerMaterial(std::vector<float>& rl)
    {
      detector_material = rl;
      detector_scatter.resize(rl.size(), 0.);
      for(unsigned int l=0;l<rl.size();++l)
      {
        detector_scatter[l] = 1.41421356237309515*0.0136*sqrt(rl[l]);
      }
    }
    
    void setVertexResolution(float vxy, float vz)
    {
      vertex_sigma_xy = vxy;
      vertex_sigma_z = vz;
    }
    
    // magnetic field in Tesla
    void setMagField(float B){Bfield = B;Bfield_inv = 1./B;}
    
    double fitTrack(SimpleTrack3D& track, std::vector<double>& chi2_hit);
    double fitTrackLine(SimpleTrack3D& track, std::vector<double>& chi2_hit);
    
    void setChi2Cut(double c){chi2_cut=c;}
    
    bool breakRecursion(const std::vector<SimpleHit3D>& hits, const HelixRange& range);
    
    float phiError(SimpleHit3D& hit, float min_k, float max_k, float min_dzdl, float max_dzdl);
    float dzdlError(SimpleHit3D& hit, float min_k, float max_k, float min_dzdl, float max_dzdl);
    
  private:
    bool using_vertex;
    float Bfield;
    float Bfield_inv;
    std::vector<float> detector_radii;
    std::vector<float> detector_radii_inv;
    std::vector<float> detector_scatter;
    std::vector<float> detector_material;
    std::vector<float> layer_xy_resolution;
    std::vector<float> layer_z_resolution;
    double vertex_sigma_xy, vertex_sigma_z;
    double chi2_cut;
    std::set<std::vector<unsigned int> > combos;
    std::set<std::vector<unsigned int> > combos_3_pass;
    std::set<std::vector<unsigned int> > combos_3_fail;
};


#endif
