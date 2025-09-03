#ifndef PARTICLEFLOWRECO_H
#define PARTICLEFLOWRECO_H

//===========================================================
/// \file ParticleFlowReco.h
/// \brief Particle flow event reconstruction
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrack;
class RawCluster;

class ParticleFlowReco : public SubsysReco
{
 public:
  ParticleFlowReco(const std::string &name = "ParticleFlowReco");

  ~ParticleFlowReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  void set_energy_match_Nsigma(float Nsigma)
  {
    _energy_match_Nsigma = Nsigma;
  }
  void set_track_map_name(std::string &name) { _track_map_name = name; }

  void set_only_crossing_zero(bool b) { _only_crossing_zero = b; }

 private:
  static int CreateNode(PHCompositeNode *topNode);

  static float calculate_dR(float, float, float, float);
  std::pair<float, float> get_expected_signature(int);

  bool _only_crossing_zero {true};

  float _energy_match_Nsigma {1.5};

  std::vector<float> _pflow_TRK_p;
  std::vector<float> _pflow_TRK_eta;
  std::vector<float> _pflow_TRK_phi;
  std::vector<float> _pflow_TRK_EMproj_phi;
  std::vector<float> _pflow_TRK_EMproj_eta;
  std::vector<float> _pflow_TRK_HADproj_phi;
  std::vector<float> _pflow_TRK_HADproj_eta;
  std::vector<SvtxTrack *> _pflow_TRK_trk;
  std::vector<std::vector<int> > _pflow_TRK_match_EM;
  std::vector<std::vector<int> > _pflow_TRK_match_HAD;

  // convention is ( EM index, dR value )
  std::vector<std::vector<std::pair<int, float> > > _pflow_TRK_addtl_match_EM;

  std::vector<float> _pflow_EM_E;
  std::vector<float> _pflow_EM_eta;
  std::vector<float> _pflow_EM_phi;
  std::vector<RawCluster *> _pflow_EM_cluster;
  std::vector<std::vector<float> > _pflow_EM_tower_eta;
  std::vector<std::vector<float> > _pflow_EM_tower_phi;
  std::vector<std::vector<int> > _pflow_EM_match_HAD;
  std::vector<std::vector<int> > _pflow_EM_match_TRK;

  std::vector<float> _pflow_HAD_E;
  std::vector<float> _pflow_HAD_eta;
  std::vector<float> _pflow_HAD_phi;
  std::vector<RawCluster *> _pflow_HAD_cluster;
  std::vector<std::vector<float> > _pflow_HAD_tower_eta;
  std::vector<std::vector<float> > _pflow_HAD_tower_phi;
  std::vector<std::vector<int> > _pflow_HAD_match_EM;
  std::vector<std::vector<int> > _pflow_HAD_match_TRK;

  std::string _track_map_name {"SvtxTrackMap"};
};

#endif  // PARTICLEFLOWRECO_H
