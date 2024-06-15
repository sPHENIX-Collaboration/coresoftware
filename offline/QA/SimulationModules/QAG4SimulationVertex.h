// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONVERTEX_H
#define QAG4SIMULATIONVERTEX_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>

#include <memory>
#include <set>
#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class TrackSeed;

class QAG4SimulationVertex : public SubsysReco
{
 public:
  QAG4SimulationVertex(const std::string &name = "QAG4SimulationVertex");

  virtual ~QAG4SimulationVertex() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void check_embed() { m_checkembed = true; }
  void embed_id_cut(const int id) { m_embed_id_cut = id; }
  void addEmbeddingID(int embeddingID);

  void setTrackgMapName(const std::string &name) { m_trackMapName = name; }
  void setVertexMapName(const std::string &name) { m_vertexMapName = name; }

 private:
  int load_nodes(PHCompositeNode *);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  unsigned int _nlayers_maps = 3;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  int m_embed_id_cut = 0;
  bool m_checkembed = false;

  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  PHG4TruthInfoContainer *m_truthInfo = nullptr;

  std::set<int> m_embeddingIDs;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_vertexMapName = "SvtxVertexMap";
};

#endif  // QAG4SIMULATIONVERTEX_H
