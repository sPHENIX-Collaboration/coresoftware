// TPC_LASERCLUSTERTRUTHMATCHER_H
#ifndef TPC_LASERCLUSTERTRUTHMATCHER_H
#define TPC_LASERCLUSTERTRUTHMATCHER_H

#include <cdbobjects/CDBTTree.h>
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
 * If wanting to run with PHGarfield distortion corrections, turn on
 * set_usePHGarfieldDistortions, and the CDB variables for PHGarfield will
 * automatically be pulled and used to propagate the clusters through the
 * fields along the poly lines before doing the matching.
 */

struct TruthRow
{
  double R{0.0};
  double dphi{0.0};
  int nstripes{0};
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
  ~LaserClusterTruthMatcher() override { delete m_cdbttree; }
  LaserClusterTruthMatcher(const LaserClusterTruthMatcher&) = delete;
  LaserClusterTruthMatcher& operator=(const LaserClusterTruthMatcher&) = delete;
  LaserClusterTruthMatcher(LaserClusterTruthMatcher&&) = delete;
  LaserClusterTruthMatcher& operator=(LaserClusterTruthMatcher&&) = delete;

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
  void setLaserClusterNodeNameOut(const std::string &name) { m_laserClusterNodeNameOut = name; }

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
  std::string m_laserClusterNodeName{"LASER_AGGREGATED_CLUSTER"};
  std::string m_laserClusterNodeNameOut{"LASER_AGGREGATED_CLUSTER_MATCHED"};
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