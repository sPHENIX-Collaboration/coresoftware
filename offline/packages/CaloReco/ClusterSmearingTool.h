// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALORECO_CLUSTERSMEARINGTOOL_H
#define CALORECO_CLUSTERSMEARINGTOOL_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <memory>  // for unique_ptr
#include <string>

class PHCompositeNode;
class RawCluster;
class RawClusterContainer;

/**
 * Applies additional energy and position smearing to (MC) clusters so that
 * the resolution matches data. The smearing model is
 *
 *   sigma_E / E   = A/E (+) B/sqrt(E) (+) C
 *   sigma_pos     = D/sqrt(E) (+) F        (applied to both eta and phi)
 *
 * where (+) denotes addition in quadrature and E is the cluster energy in GeV.
 * The parameters are read from the CDB (calib name CEMC_ClusterSmearing by
 * default) or from a file given with set_directURL().
 *
 * By default the clusters in the input node are modified in place. Calling
 * set_outputNodeName() instead writes smeared copies into a new node and leaves
 * the input node untouched.
 */
class ClusterSmearingTool : public SubsysReco
{
 public:
  explicit ClusterSmearingTool(const std::string &name = "ClusterSmearingTool");

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  void setCalibName(const std::string &name)
  {
    m_calibName = name;
    return;
  }

  void set_directURL(const std::string &url)
  {
    m_giveDirectURL = true;
    m_directURL = url;
  }

  void set_inputNodeName(const std::string &name)
  {
    m_inputNodeName = name;
    return;
  }

  // write smeared clusters to a new node instead of editing the input in place
  void set_outputNodeName(const std::string &name)
  {
    m_outputNodeName = name;
    m_writeNewNode = true;
    return;
  }

  void set_seed(unsigned int seed)
  {
    m_seed = seed;
    m_overrideSeed = true;
  }

  void set_doEnergySmearing(bool doit) { m_doEnergySmearing = doit; }
  void set_doPositionSmearing(bool doit) { m_doPositionSmearing = doit; }

  // eta smearing is done relative to the event vertex (GlobalVertexMap) when available
  void set_useVertex(bool use) { m_useVertex = use; }

  void set_doAbortMissingCalib(bool doAbort = true)
  {
    m_doAbortMissingCalib = doAbort;
    return;
  }

  float get_A() const { return m_parA; }
  float get_B() const { return m_parB; }
  float get_C() const { return m_parC; }
  float get_D() const { return m_parD; }
  float get_F() const { return m_parF; }

 private:
  void LoadCalib();
  void SmearCluster(RawCluster *cluster, float vtx_z);

  std::string m_detector{"CEMC"};
  std::string m_calibName{"CEMC_ClusterSmearing"};
  std::string m_inputNodeName{"CLUSTERINFO_CEMC"};
  std::string m_outputNodeName;
  bool m_writeNewNode{false};

  bool m_giveDirectURL{false};
  std::string m_directURL;

  bool m_doAbortMissingCalib{true};
  bool m_doSmearing{true};
  bool m_doEnergySmearing{true};
  bool m_doPositionSmearing{true};
  bool m_useVertex{true};

  unsigned int m_seed{0};
  bool m_overrideSeed{false};

  // smearing parameters, see class description
  float m_parA{0};
  float m_parB{0};
  float m_parC{0};
  float m_parD{0};
  float m_parF{0};

  class Deleter
  {
   public:
    void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
  };
  std::unique_ptr<gsl_rng, Deleter> m_rng;
};

#endif  // CALORECO_CLUSTERSMEARINGTOOL_H
