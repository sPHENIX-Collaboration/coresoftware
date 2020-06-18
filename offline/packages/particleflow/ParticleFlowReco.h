#ifndef PARTICLEFLOWRECO_H
#define PARTICLEFLOWRECO_H

//===========================================================
/// \file ParticleFlowReco.h
/// \brief Particle flow event reconstruction
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;

class ParticleFlowReco : public SubsysReco
{
 public:

  ParticleFlowReco(const std::string &name = "ParticleFlowReco");

  virtual ~ParticleFlowReco();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 private:

  int CreateNode(PHCompositeNode *topNode);

  float calculate_dR( float, float, float, float );
  std::pair<float, float> get_expected_signature( int );

  std::vector<float> _pflow_TRK_p;
  std::vector<float> _pflow_TRK_eta;
  std::vector<float> _pflow_TRK_phi;
  std::vector< std::vector<int> > _pflow_TRK_match_EM;
  std::vector< std::vector<int> > _pflow_TRK_match_HAD;

  // convention is ( EM index, dR value )
  std::vector< std::vector< std::pair<int,float> > > _pflow_TRK_addtl_match_EM;

  std::vector<float> _pflow_EM_E;
  std::vector<float> _pflow_EM_eta;
  std::vector<float> _pflow_EM_phi;
  std::vector< std::vector<float> > _pflow_EM_tower_eta;
  std::vector< std::vector<float> > _pflow_EM_tower_phi;
  std::vector< std::vector<int> > _pflow_EM_match_HAD;
  std::vector< std::vector<int> > _pflow_EM_match_TRK;

  std::vector<float> _pflow_HAD_E;
  std::vector<float> _pflow_HAD_eta;
  std::vector<float> _pflow_HAD_phi;
  std::vector< std::vector<float> > _pflow_HAD_tower_eta;
  std::vector< std::vector<float> > _pflow_HAD_tower_phi;
  std::vector< std::vector<int> > _pflow_HAD_match_EM;
  std::vector< std::vector<int> > _pflow_HAD_match_TRK;


};

#endif // PARTICLEFLOWRECO_H
