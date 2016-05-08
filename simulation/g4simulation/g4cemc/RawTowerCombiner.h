#ifndef RawTowerCombiner_H__
#define RawTowerCombiner_H__

#include <fun4all/SubsysReco.h>
#include <string>

#include <phool/PHTimeServer.h>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

//! RawTowerCombiner module that joints multiple RawTower together to form a single readout in a separate node
//! Use this class to simulate ganged readout or trigger channels that combined multiple towers
//! Should be called after RawTowerBuilder but before RawTowerDigitizer
class RawTowerCombiner : public SubsysReco
{

public:
  RawTowerCombiner(const std::string& name = "RawTowerCombiner");

  virtual
  ~RawTowerCombiner()
  {
  }

  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);
  void
  Detector(const std::string &d)
  {
    detector = d;
  }

  std::string
  get_tower_node_prefix() const
  {
    return _tower_node_prefix;
  }

  void
  set_tower_node_prefix(std::string simTowerNodePrefix)
  {
    _tower_node_prefix = simTowerNodePrefix;
  }

  unsigned int
  get_combine_eta() const
  {
    return _n_combine_eta;
  }

  void
  set_combine_eta(unsigned int combineEta)
  {
    _n_combine_eta = combineEta;
  }

  unsigned int
  get_combine_phi() const
  {
    return _n_combine_phi;
  }

  void
  set_combine_phi(unsigned int combinePhi)
  {
    _n_combine_phi = combinePhi;
  }

  std::string
  get_output_node_suffix() const
  {
    return _output_node_suffix;
  }

  void
  set_output_node_suffix(std::string outputNodeSuffix)
  {
    _output_node_suffix = outputNodeSuffix;
  }
protected:

  std::string _tower_node_prefix;

  //! if empty, suffix will be automatically assigned
  std::string _output_node_suffix;

  unsigned int _n_combine_eta;
  unsigned int _n_combine_phi;

  //! get the new ieta from the old
  inline int get_output_bin_eta(int input_bin) const {return input_bin/_n_combine_eta;}
  //! get the new iphi from the old
  inline int get_output_bin_phi(int input_bin) const {return input_bin/_n_combine_phi;}


  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _intput_towers;
  RawTowerContainer* _output_towers;

  std::string detector;
  std::string TowerNodeName;
  std::string TowerGeomNodeName;

};

#endif /* RawTowerCombiner_H__ */
