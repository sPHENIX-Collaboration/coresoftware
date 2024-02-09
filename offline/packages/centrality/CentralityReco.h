#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>
#include <array>
#include <limits>
#include <string>  // for string, allocator

// Forward declarations
class CentralityInfo;
class PHCompositeNode;
class MbdOut;

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  explicit CentralityReco(const std::string &name = "CentralityReco");

  //! destructor
  ~CentralityReco() override = default;

  //! full initialization
  int InitRun(PHCompositeNode *) override;
  void CreateNodes(PHCompositeNode *);
  int FillCentralityInfo();
  int FillVars();
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *) override;
  int ResetEvent(PHCompositeNode *) override;

  // Interface with CDB
  int Download_centralityDivisions(const std::string &dbfile);
  int Download_centralityScale(const std::string &dbfile);

 private:
  std::string _dbfilename;

  const int NDIVS{100};

  MbdOut *_mbd_out{nullptr};

  CentralityInfo *_central{nullptr};

  unsigned int _key{std::numeric_limits<unsigned int>::max()};

  float _mbd_charge_sum{std::numeric_limits<float>::quiet_NaN()};
  float _mbd_charge_sum_n{std::numeric_limits<float>::quiet_NaN()};
  float _mbd_charge_sum_s{std::numeric_limits<float>::quiet_NaN()};

  double _centrality_scale{std::numeric_limits<double>::quiet_NaN()};
  std::array<float, 100> _centrality_map{};

};

#endif
