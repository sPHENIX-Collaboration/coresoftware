// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file		PHG4TrackFastSimEval.h
 *  \brief		Evaluation module for PHG4TrackFastSim output
 *  \details	input: PHG4TruthInfoContainer, SvtxTrackMap with SvtxTrack_FastSim inside
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H
#define G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H

#include <fun4all/SubsysReco.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

// Forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class TTree;
class TH2;

// Brief: basic ntuple and histogram creation for sim evaluation
class PHG4TrackFastSimEval : public SubsysReco
{
 public:
  // Default constructor
  PHG4TrackFastSimEval(const std::string& name = "PHG4TrackFastSimEval",
                       const std::string& filename = "g4eval.root",
                       const std::string& trackmapname = "SvtxTrackMap");

  // Initialization, called for initialization
  int Init(PHCompositeNode*) override;

  // Initialization, called for initialization
  int InitRun(PHCompositeNode*) override;

  // Process Event, called for each event
  int process_event(PHCompositeNode*) override;

  // End, write and close files
  int End(PHCompositeNode*) override;

  // Change output filename
  void set_filename(const std::string& file)
  {
    m_OutFileName = file;
  }

  // set the name of the node with the trackmap
  void set_trackmapname(const std::string& name)
  {
    m_TrackMapName = name;
  }

  // User modules
  void reset_variables();

  void AddProjection(const std::string& name);

 private:
  void fill_track_tree(PHCompositeNode*);
  void fill_vertex_tree(PHCompositeNode*);

  // Get all the nodes
  int GetNodes(PHCompositeNode*);

  // Node pointers
  PHG4TruthInfoContainer* m_TruthInfoContainer{nullptr};
  SvtxTrackMap* m_TrackMap{nullptr};
  SvtxVertexMap* m_VertexMap{nullptr};

  // TTrees
  TTree* m_TracksEvalTree{nullptr};
  TTree* m_VertexEvalTree{nullptr};

  // Histos
  TH2* m_H2D_DeltaMomVsTruthMom{nullptr};
  TH2* m_H2D_DeltaMomVsTruthEta{nullptr};

  // Event counter
  int m_EventCounter{0};

  // TTree variables
  int m_TTree_Event{-9999};
  //-- truth
  int m_TTree_gTrackID{-9999};
  int m_TTree_gFlavor{-9999};
  float m_TTree_gpx{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gpy{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gpz{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gvx{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gvy{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gvz{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_gvt{std::numeric_limits<float>::quiet_NaN()};

  //-- reco
  int m_TTree_TrackID{-9999};
  int m_TTree_Charge{-9999};
  int m_TTree_nHits{-9999};
  float m_TTree_px{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_py{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_pz{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_pcax{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_pcay{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_pcaz{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_dca2d{std::numeric_limits<float>::quiet_NaN()};

  std::map<int, int> m_TTree_HitContainerID_nHits_map;

  // vertex
  float m_TTree_vx{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_vy{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_vz{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_DeltaVx{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_DeltaVy{std::numeric_limits<float>::quiet_NaN()};
  float m_TTree_DeltaVz{std::numeric_limits<float>::quiet_NaN()};
  int m_TTree_nTracks{-9999};
  int m_TTree_nFromTruth{-9999};

  // output filename
  std::string m_OutFileName;

  // name of SvtxTrackMap collection
  std::string m_TrackMapName;

  // names and index of projections
  std::map<std::string, unsigned int> m_ProjectionNameMap;
  // projections to cylinders and planes
  std::vector<std::vector<float>> m_TTree_proj_vec;
  std::vector<std::vector<float>> m_TTree_proj_p_vec;
  // hits on reference cylinders and planes
  std::vector<std::vector<float>> m_TTree_ref_vec;
  std::vector<std::vector<float>> m_TTree_ref_p_vec;
};

#endif  //* G4TRACKFASTSIM_PHG4TRACKFASTSIMEVAL_H *//
