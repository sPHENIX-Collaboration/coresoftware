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

  //! prefix to the tower node
  std::string
  get_tower_node_prefix() const
  {
    return _tower_node_prefix;
  }

  //! prefix to the tower node
  void
  set_tower_node_prefix(const std::string & simTowerNodePrefix)
  {
    _tower_node_prefix = simTowerNodePrefix;
  }

  //! number of eta towers to be merged into a new tower
  unsigned int
  get_combine_eta() const
  {
    return _n_combine_eta;
  }

  //! number of eta towers to be merged into a new tower
  void
  set_combine_eta(unsigned int combineEta)
  {
    _n_combine_eta = combineEta;
  }

  //! number of eta towers to be merged into a new tower
  unsigned int
  get_combine_phi() const
  {
    return _n_combine_phi;
  }

  //! number of eta towers to be merged into a new tower
  void
  set_combine_phi(unsigned int combinePhi)
  {
    _n_combine_phi = combinePhi;
  }

protected:

  //! prefix to the tower node
  std::string _tower_node_prefix;

  //! number of eta towers to be merged into a new tower
  unsigned int _n_combine_eta;
  //! number of phi towers to be merged into a new tower
  unsigned int _n_combine_phi;

  //! get the new ieta from the old
  inline int get_output_bin_eta(int input_bin) const {return input_bin/_n_combine_eta;}
  //! get the new iphi from the old
  inline int get_output_bin_phi(int input_bin) const {return input_bin/_n_combine_phi;}


  void
  CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer* _towers;

  std::string detector;

};

#endif /* RawTowerCombiner_H__ */
