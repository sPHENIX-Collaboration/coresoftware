#ifndef TPC_MODULETRACKDISPLAY_H
#define TPC_MODULETRACKDISPLAY_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>

class PHCompositeNode;
class TFile;
class Tpc_ModuleTrackContainer;
class TrkrHitSetContainer;
class IdealPadMap;

class Tpc_ModuleTrackDisplay : public SubsysReco
{
 public:
  Tpc_ModuleTrackDisplay(const std::string& name = "Tpc_ModuleTrackDisplay",
                       const std::string& outfilename = "tpc_moduletrack_display.root",
                       const std::string& trackNodeName = "TPC_MODULETRACKS",
                       unsigned int maxEventDisplays = 5);
  ~Tpc_ModuleTrackDisplay() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  struct HitPoint
  {
    HitPoint();

    bool ok;
    TrkrDefs::hitsetkey hitsetkey;
    TrkrDefs::hitkey hitkey;

    int side;
    unsigned int sector;
    unsigned int layer;
    unsigned int pad;
    unsigned int tbin;
    unsigned short adc;

    // Only display/fitter coordinates.  No global/local x/y are used.
    double radius;
    double phi;
  };

  // Display-only fit mode. Use Tpc_FittingTools::FIT_SAGITTA for field-on curvature
  // or Tpc_FittingTools::FIT_LINEAR for straight-line fits. Default is sagitta.
  void setFitMode(int mode) { m_fitMode = mode; }
  void setFitWithSagitta(bool v) { m_fitMode = v ? 1 : 0; }

  void setFitWeightPower(double v) { m_fitWeightPower = v; }
  void setFitWeightFloorFrac(double v) { m_fitWeightFloorFrac = v; }

 private:
  bool get_nodes(PHCompositeNode* topNode);
  HitPoint make_hit_point(TrkrDefs::hitsetkey hsk,
                          TrkrDefs::hitkey hk) const;

  double ideal_radius(unsigned int layer) const;
  double ideal_phi(int side,
                   unsigned int sector,
                   unsigned int layer,
                   unsigned int pad) const;

  std::string m_outfilename;
  std::string m_trackNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;

  TFile* m_outfile;
  Tpc_ModuleTrackContainer* m_tracks;
  TrkrHitSetContainer* m_hits;
  IdealPadMap* m_idealPadMap;

  int m_fitMode;
  double m_fitWeightPower;
  double m_fitWeightFloorFrac;
};

#endif
