#ifndef TRACKBASEHISTORIC_TRACKSEEDHELPER_H
#define TRACKBASEHISTORIC_TRACKSEEDHELPER_H

class TrackSeed;
class TrkrClusterContainer;
class ActsGeometry;

class TrackSeedHelper
{

  public:

  float get_px( TrackSeed*, TrkrClusterContainer*, ActsGeometry* );
  float get_py( TrackSeed*, TrkrClusterContainer*, ActsGeometry* );
  float get_phi( TrackSeed*, TrkrClusterContainer*, ActsGeometry* );
  float get_phi(const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions) const override;

  float get_phi_fastsim( TrackSeed*, TrkrClusterContainer*, ActsGeometry* );

  void circleFitByTaubin(TrkrClusterContainer* clusters,
    ActsGeometry* tGeometry,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58) override;

  void lineFit(TrackSeed*, TrkrClusterContainer* clusters,
    ActsGeometry* tGeometry,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58) override;

  void circleFitByTaubin(TrackSeed*, const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58) override;

  void lineFit(TrackSeed*, const std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58) override;

  float get_x(TrackSeed*) const override;
  float get_y(TrackSeed*) const override;
  float get_z(TrackSeed*) const override;


  protected:
  std::pair<float, float> findRoot(TrackSeed*) const;

};



#endif
