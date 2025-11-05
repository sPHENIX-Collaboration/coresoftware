#ifndef CALOVTXRECO_H
#define CALOVTXRECO_H
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <globalvertex/CaloVertexMapv1.h>
class PHCompositeNode;
class CaloVtxReco : public SubsysReco
{
 public:

  CaloVtxReco(const std::string &name = "CaloVtxReco", const std::string &jetnodename = "zzjets06", const int debug = 0, const bool use_z_energy_dep = false);

  virtual ~CaloVtxReco();

  int createNodes(PHCompositeNode *topNode);

  float new_eta(int channel, TowerInfoContainer* towers, RawTowerGeomContainer* geom, RawTowerDefs::CalorimeterId caloID, float testz);
  
  int InitRun(PHCompositeNode *topNode) override;

  int calo_tower_algorithm(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  void Print(const std::string &what = "ALL") const override;

  float get_jet_threshold() { return _jet_threshold; }

  void set_jet_threshold(float new_thresh) { _jet_threshold = new_thresh; }

  float get_calib_factor() { return _calib_factor; }

  void set_calib_factor(float new_calib) { _calib_factor = new_calib; }

  float get_energy_cut() { return _energy_cut; }

  void set_energy_cut(float new_energy) { _energy_cut = new_energy; }
 
 private:
  std::string _name;
  std::string _jetnodename;
  int _debug;
  bool _use_z_energy_dep;
  float _jet_threshold{15};
  float _zvtx{-9999};
  float _calib_factor{1.406};
  float _energy_cut{0.1};
  CaloVertexMap* _calovtxmap{};
};

#endif // CALOVTXRECO
