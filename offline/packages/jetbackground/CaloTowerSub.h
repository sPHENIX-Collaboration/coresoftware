#ifndef JETBACKGROUND_CALOTOWERSUB_H
#define JETBACKGROUND_CALOTOWERSUB_H

#include <fun4all/SubsysReco.h>

#include <jetbase/Jet.h>

#include <string>
#include <vector>

#include <globalvertex/GlobalVertex.h>

class PHCompositeNode;

class JetAlgo;
class JetInput;

class CaloTowerSub : public SubsysReco
{
 
 public:
 
  CaloTowerSub( const std::string &name = "CaloTowerSub" );
  ~CaloTowerSub() override;

  int InitRun( PHCompositeNode *topNode ) override;
  int process_event( PHCompositeNode *topNode ) override;

  enum CALOSUB_PARAM 
  { 
    SEED_MIN_D = 0,
    SEED_MIN_LEAD_COMP_ET = 1,
    SEED_MIN_RAW_PT = 2,
    SEED_MIN_SUB_PT = 3,
    EXCLUSION_DR = 4,
    FLOW_MODE = 5,
    REWEIGHT = 6
  };

  void setParam(
    const CALOSUB_PARAM opt,
    float value
  );

  void setParam(
    const char * opt,
    float value
  );

  void addCaloInput( JetInput * input ) 
  { 
    m_inputs.push_back(input); 
  }

  void setSeedAlgo( JetAlgo * algo, const std::string & output)
  {
    m_seed_algo = algo;
    m_raw_seed_output = output;
  }



 private:

  int CreateNodes(PHCompositeNode *topNode);
  void FillNodes(PHCompositeNode *topNode);

  std::vector<JetInput *> m_inputs {};
  JetAlgo *m_seed_algo {nullptr};
  std::string m_raw_seed_output {""};
  int recoSeeds( PHCompositeNode *topNode , const int iter );

  float m_zvrtx {0.0};
  GlobalVertex::VTXTYPE m_vertex_type = GlobalVertex::UNDEFINED;
  int getZvertex( PHCompositeNode *topNode );
  
  std::vector< std::vector < std::vector< int > > > _CALO_TOWER_MASK {}; // layer, ieta, iphi
  std::vector< std::vector < std::vector< float > > > _CALO_E {}; // layer, ieta, iphi
  std::vector< std::vector < std::vector< int > > > _CALO_ISBAD {}; // layer, ieta, iphi
  std::vector< std::vector < std::vector< std::pair< float, float > > > > _CALO_POS {}; // layer, ieta, iphi
  
  std::vector< std::vector < std::vector< float > > > _CALO_WEIGHTS {}; // layer, ieta, iphi
  std::vector< std::vector < std::vector< float > > > _CALO_FLOW_WEIGHTS {}; // layer, ieta, iphi

  int genSeedMask( PHCompositeNode *topNode );
  int   m_seed_type{0};
  float m_seed_jet_D{4.0};
  float m_seed_max_const{3.0};
  float m_seed_jet_pt[2] = {5.0, 7.0};
  float m_seed_exclusion_dR{0.4};

  std::vector < float > _seed_eta;
  std::vector < float > _seed_phi;

  Jet::PROPERTY _index_SeedD{};
  Jet::PROPERTY _index_SeedItr{};

  int _HCAL_NETA {-1};
  int _HCAL_NPHI {-1};
  int _N_LAYERS  {-1};

  

  // std::vector< std::vector < std::vector< int > > > _CALO_TOWER_MASK {}; // layer, ieta, iphi
  // std::vector< std::vector < std::vector< float > > > _CALO_E {}; // layer, ieta, iphi
  // std::vector< std::vector < std::vector< int > > > _CALO_ISBAD {}; // layer, ieta, iphi
  // std::vector< std::vector < std::vector< std::pair< float, float > > > > _CALO_POS {}; // layer, ieta, iphi
  
  // std::vector< std::vector < std::vector< float > > > _CALO_WEIGHTS {}; // layer, ieta, iphi
  // std::vector< std::vector < std::vector< float > > > _CALO_FLOW_WEIGHTS {}; // layer, ieta, iphi
  
  
  int fillCaloArrays( PHCompositeNode *topNode );  
  int getCaloNodes( PHCompositeNode *topNode );
  int calcFlow( PHCompositeNode *topNode );
  int calcUE( PHCompositeNode *topNode );

  int _do_flow{0};
  float _v2{0};
  float _Psi2{0};
  std::vector<std::vector<float> > _UE;
  int _nStrips{0};
  int _nTowers{0};



  bool getCaloStatus( const int layer, const int ieta, const int iphi );
  float getCaloEta( const int layer, const int ieta, const int iphi );
  float getCaloPhi( const int layer, const int ieta, const int iphi );
  
  static float deltaR( const float eta1, const float phi1, const float eta2, const float phi2 );

  std::vector<std::vector<float> > _EMCAL_E;
  std::vector<std::vector<float> > _IHCAL_E;
  std::vector<std::vector<float> > _OHCAL_E;

  std::vector<std::vector<int> > _EMCAL_ISBAD;
  std::vector<std::vector<int> > _IHCAL_ISBAD;
  std::vector<std::vector<int> > _OHCAL_ISBAD;

  // 1-D energies vs. phi (integrated over eta strips with complete
  // phi coverage, and all layers)
  std::vector<float> _FULLCALOFLOW_PHI_E;
  std::vector<float> _FULLCALOFLOW_PHI_VAL;

  bool _do_reweight{true}; // flag to indicate if reweighting is used
  std::vector<float> _EMCAL_PHI_WEIGHTS;
  std::vector<float> _IHCAL_PHI_WEIGHTS;
  std::vector<float> _OHCAL_PHI_WEIGHTS;

  std::string _backgroundName{"TestTowerBackground"};

  

  bool m_use_towerinfo{false};
  bool _is_flow_failure{false};
  bool _reweight_failed{false};

  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;
};

#endif
