// TODO: put this in the right coresoftware package (probably tpc/, next to
// LaserPadAggregator / LaserAggregatedClusterizer) and fix the include guard
// / header path below to match.
#ifndef TPC_LASERCLUSTERTRUTHMATCHER_H
#define TPC_LASERCLUSTERTRUTHMATCHER_H

#include <tpc/LaserClusterHelper.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <vector>


class PHCompositeNode;
class LaserClusterContainer;
class CDBTTree;

/*
 * Matches already-aggregated-and-clustered laser clusters (LaserClusterv5,
 * produced upstream by LaserAggregatedClusterizer) to the CDB truth laser
 * pattern, writing only a truthIndex (or -999 if unmatched) back onto each
 * cluster in place.
 *
 * This is deliberately a separate module from LaserAggregatedClusterizer:
 * clustering is per-block/parallelized, matching needs the FULL set of
 * clusters for a side at once (row-finding and stripe-finding are whole-side
 * KDE operations), and matching needs to be retuned/rerun far more often
 * than clustering itself. Because distortions differ run to run, matching is
 * always redone from scratch here -- nothing is cached across runs.
 *
 * This module never loads its own TpcDistortionCorrectionContainer -- it
 * expects one to already be on the node tree (loaded by the macro / an
 * upstream module) and will abort the run if it isn't there.
 */

struct TruthRow
{
  double R;
  double dphi;
  int nstripes;
};

struct TruthRowPattern
{
  std::vector<double> stripePhi;   // coordinate phi, ascending
  std::vector<int> stripeIPhi;          // stripeIPhi[slot] = CDB iphi digit for that slot
};

class LaserClusterTruthMatcher : public SubsysReco
{
 public:
  explicit LaserClusterTruthMatcher(const std::string &name = "LaserClusterTruthMatcher");
  ~LaserClusterTruthMatcher() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //! path to the CDB truth-pattern file (CMStripePattern_*.root)
  //  TODO: this should probably come from the CDB machinery (run-dependent
  //  calibration lookup) rather than a hardcoded path, once this is wired
  //  into the real calibration framework.
  void setTruthFile(const std::string &file) { m_truthFile = file; }
  
  //! node name of the LaserClusterContainer written by LaserAggregatedClusterizer
  void setLaserClusterNodeName(const std::string &name) { m_laserClusterNodeName = name; }

  void set_useGlobal(bool use) { m_useGlobal = use; }

  void set_usePHGarfieldDistortions(bool use) { m_usePHGarfieldDistortions = use; }
  void set_garfield_cmvoltage(double use) { m_garfield_cmvoltage = use; }
  void set_garfield_zerofield(bool use) { m_garfield_zerofield = use; }
  void set_garfield_keffside0(double use) { m_garfield_keffside0 = use; }
  void set_garfield_keffside1(double use) { m_garfield_keffside1 = use; }
  void set_garfield_stepns(double use) { m_garfield_stepns = use; }

  void set_QABase(const std::string &name) { m_QABase = name; }


 private:
  int getNodes(PHCompositeNode *topNode);

  std::string m_truthFile{"CMStripePattern_full.root"};
  std::string m_laserClusterNodeName{"LASER_CLUSTER"};  // TODO: confirm actual node name
  std::string m_QABase{""};

  LaserClusterContainer *m_laserClusterContainer{nullptr};
  CDBTTree *m_cdbttree{nullptr};

  std::map<int, TruthRowPattern> m_truthRowPatterns[2];
  std::vector<TruthRow> m_truthRows[2];

  LaserClusterHelper m_laserClusterHelper;
  bool m_useGlobal{true};
  bool m_usePHGarfieldDistortions{false};
  double m_garfield_cmvoltage{380.0};
  bool m_garfield_zerofield{false};
  double m_garfield_keffside0{1.0};
  double m_garfield_keffside1{1.0};
  double m_garfield_stepns{50.0};
  
};

#endif