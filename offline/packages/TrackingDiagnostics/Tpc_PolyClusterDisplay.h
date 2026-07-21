#ifndef TPC_POLYCLUSTERDISPLAY_H
#define TPC_POLYCLUSTERDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Tpc_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;
class TFile;
class Tpc_PolyClusterContainer;

class Tpc_PolyClusterDisplay : public SubsysReco
{
 public:
  Tpc_PolyClusterDisplay(const std::string& name = "Tpc_PolyClusterDisplay",
                        const std::string& outfilename = "tpc_polycluster_display.root",
                        const std::string& clusterNodeName = "TPC_POLYCLUSTERS",
                        unsigned int maxEventDisplays = 5);
  ~Tpc_PolyClusterDisplay() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setTpc_PolyTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setTpc_PolyTrackVertexNodeName(const std::string& n) { m_finalTrackVertexNodeName = n; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setZRange(double zmin, double zmax) { m_zmin = zmin; m_zmax = zmax; }
  void setTrackVertexZRange(double zmin, double zmax) { m_trackVertexZMin = zmin; m_trackVertexZMax = zmax; }
  void setXYRange(double xymax) { m_xymax = xymax; }
  void setUseStraightLineTracks(bool v) { m_useStraightLineTracks = v; }

 private:
  bool get_nodes(PHCompositeNode* topNode);

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_finalTrackNodeName;
  std::string m_finalTrackVertexNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;
  double m_zmin;
  double m_zmax;
  double m_trackVertexZMin;
  double m_trackVertexZMax;
  double m_xymax;
  double m_magneticFieldTesla;
  bool m_useStraightLineTracks;

  TFile* m_outfile;
  Tpc_PolyClusterContainer* m_clusters;
  Tpc_PolyTrackContainer* m_finalTracks;
  Tpc_PolyTrackVertexContainer* m_finalTrackVertices;
};

#endif
