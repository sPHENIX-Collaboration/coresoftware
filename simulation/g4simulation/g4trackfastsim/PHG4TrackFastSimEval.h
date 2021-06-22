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

#include <map>
#include <string>
#include <vector>

//Forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class TTree;
class TH2D;

//Brief: basic ntuple and histogram creation for sim evaluation
class PHG4TrackFastSimEval : public SubsysReco
{
 public:
  //Default constructor
  PHG4TrackFastSimEval(const std::string& name = "PHG4TrackFastSimEval",
                       const std::string& filename = "g4eval.root",
                       const std::string& trackmapname = "SvtxTrackMap");

  //Initialization, called for initialization
  int Init(PHCompositeNode*) override;

  //Initialization, called for initialization
  int InitRun(PHCompositeNode*) override;

  //Process Event, called for each event
  int process_event(PHCompositeNode*) override;

  //End, write and close files
  int End(PHCompositeNode*) override;

  //Change output filename
  void set_filename(const std::string& file)
  {
    m_OutFileName = file;
  }

  // set the name of the node with the trackmap
  void set_trackmapname(const std::string& name)
  {
    m_TrackMapName = name;
  }

  //User modules
  void reset_variables();

  void AddProjection(const std::string& name);

 private:
  void fill_track_tree(PHCompositeNode*);
  void fill_vertex_tree(PHCompositeNode*);

  //Get all the nodes
  int GetNodes(PHCompositeNode*);

  //Node pointers
  PHG4TruthInfoContainer* m_TruthInfoContainer;
  SvtxTrackMap* m_TrackMap;
  SvtxVertexMap* m_VertexMap;

  //TTrees
  TTree* m_TracksEvalTree;
  TTree* m_VertexEvalTree;

  //Histos
  TH2D* m_H2D_DeltaMomVsTruthMom;
  TH2D* m_H2D_DeltaMomVsTruthEta;

  //Event counter
  int m_EventCounter;

  // TTree variables
  int m_TTree_Event;
  //-- truth
  int m_TTree_gTrackID;
  int m_TTree_gFlavor;
  float m_TTree_gpx;
  float m_TTree_gpy;
  float m_TTree_gpz;
  float m_TTree_gvx;
  float m_TTree_gvy;
  float m_TTree_gvz;
  float m_TTree_gvt;

  //-- reco
  int m_TTree_TrackID;
  int m_TTree_Charge;
  int m_TTree_nHits;
  float m_TTree_px;
  float m_TTree_py;
  float m_TTree_pz;
  float m_TTree_pcax;
  float m_TTree_pcay;
  float m_TTree_pcaz;
  float m_TTree_dca2d;

  std::map<int, int> m_TTree_HitContainerID_nHits_map;

  //vertex
  float m_TTree_vx;
  float m_TTree_vy;
  float m_TTree_vz;
  float m_TTree_DeltaVx;
  float m_TTree_DeltaVy;
  float m_TTree_DeltaVz;
  int m_TTree_nTracks;
  int m_TTree_nFromTruth;

  //output filename
  std::string m_OutFileName;

  //name of SvtxTrackMap collection
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
