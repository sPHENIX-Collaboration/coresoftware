#ifndef TRACKBASEHISTORIC_TRACKSEEDHELPER_H
#define TRACKBASEHISTORIC_TRACKSEEDHELPER_H

/*!
 * \file TrackSeedHelper.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@lanl.gov>
 */

class TrackSeed;

class TrackSeedHelper
{

  public:

  using position_map_t = std::map<TrkrDefs::cluskey, Acts::Vector3>;

  static float get_phi(TrackSeed const*, const position_map_t&);

  static void circleFitByTaubin(
    TrackSeed*, const position_map_t& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58);

  static void lineFit(
    TrackSeed*, const position_map_t& positions,
    uint8_t startLayer = 0,
    uint8_t endLayer = 58);

  static float get_x(TrackSeed const*);
  static float get_y(TrackSeed const*);
  static float get_z(TrackSeed const*);

  protected:
  static std::pair<float, float> findRoot(TrackSeed const*);

};



#endif
