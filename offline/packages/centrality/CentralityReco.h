#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>

<<<<<<< HEAD
#include <limits>

// Forward declarations
class CentralityInfo;
class Fun4AllHistoManager;
class PHCompositeNode;
class CDBInterface;
class CDBTTree;
class recoConsts;
class MbdOutV1;

=======
#include <array>
#include <limits>
#include <string>  // for string, allocator

// Forward declarations
class CentralityInfo;
class PHCompositeNode;
class MbdOut;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  explicit CentralityReco(const std::string &name = "CentralityReco");

  //! destructor
<<<<<<< HEAD
  virtual ~CentralityReco();

  //! full initialization
  int Init(PHCompositeNode *) override;
=======
  ~CentralityReco() override = default;

  //! full initialization
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  int InitRun(PHCompositeNode *) override;
  void CreateNodes(PHCompositeNode *);
  int FillCentralityInfo();
  int FillVars();
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *) override;
<<<<<<< HEAD
  
  //! end of run method
  int End(PHCompositeNode *) override;
=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  int ResetEvent(PHCompositeNode *) override;

  // Interface with CDB
<<<<<<< HEAD
  int Download_centralityDivisions(const std::string& dbfile);
  int Download_centralityScale(const std::string& dbfile);

 private:

  CDBInterface *_cdb {nullptr};
  recoConsts   *_rc {nullptr};
  std::string  _dbfilename;

  const int NDIVS = 18;

  MbdOutV1 *_mbd_out = nullptr;

  CentralityInfo *_central = nullptr;

  unsigned int _key = std::numeric_limits<unsigned int>::max();

  float _mbd_charge_sum = std::numeric_limits<float>::quiet_NaN();
  float _mbd_charge_sum_n = std::numeric_limits<float>::quiet_NaN();
  float _mbd_charge_sum_s = std::numeric_limits<float>::quiet_NaN();

  double _centrality_scale = std::numeric_limits<double>::quiet_NaN();
  float _centrality_map[20]{};
=======
  int Download_centralityDivisions(const std::string &dbfile);
  int Download_centralityScale(const std::string &dbfile);

 private:
  std::string _dbfilename;

  const int NDIVS{18};

  MbdOut *_mbd_out{nullptr};

  CentralityInfo *_central{nullptr};

  unsigned int _key{std::numeric_limits<unsigned int>::max()};

  float _mbd_charge_sum{std::numeric_limits<float>::quiet_NaN()};
  float _mbd_charge_sum_n{std::numeric_limits<float>::quiet_NaN()};
  float _mbd_charge_sum_s{std::numeric_limits<float>::quiet_NaN()};

  double _centrality_scale{std::numeric_limits<double>::quiet_NaN()};
  std::array<float, 20> _centrality_map{};
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
};

#endif
