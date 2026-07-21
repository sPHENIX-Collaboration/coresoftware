#pragma once

#include "Tpc_FittingTools.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class Tpc_PolyTrackContainer;
class IdealPadMap;
class PHCompositeNode;
class Tpc_PolyCluster;
class Tpc_PolyClusterContainer;

class Tpc_PolyTrackReco : public SubsysReco
{
 public:
  enum class FitMode
  {
    Helix,
    Line3D
  };

  explicit Tpc_PolyTrackReco(const std::string& name = "Tpc_PolyTrackReco");
  ~Tpc_PolyTrackReco() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }
  void setFitMode(FitMode mode) { m_fitMode = mode; }
  void setUseLine3DFit(bool v) { m_fitMode = v ? FitMode::Line3D : FitMode::Helix; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  double calc_dedx(const std::vector<const Tpc_PolyCluster*>& clusters,
                   const Tpc_FittingTools::FitResult& fit,
                   bool fit_ok) const;
  void fillTpc_PolyTrack(unsigned int source_assembled_track_id,
                         const std::vector<const Tpc_PolyCluster*>& clusters,
                      const Tpc_FittingTools::FitResult& fit,
                      bool fit_ok);

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  Tpc_PolyClusterContainer* m_clusters {nullptr};
  Tpc_PolyTrackContainer* m_polyTracks {nullptr};
  IdealPadMap* m_idealPadMap {nullptr};
  unsigned int m_event {0};
  double m_magneticFieldTesla {1.4};
  FitMode m_fitMode {FitMode::Helix};
};
