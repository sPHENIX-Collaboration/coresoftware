// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_CDBINTERFACE_H
#define FFAMODULES_CDBINTERFACE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>  // for uint64_t
#include <set>
#include <string>
#include <tuple>  // for tuple

class PHCompositeNode;
class SphenixClient;

class CDBInterface : public SubsysReco
{
 public:
  static CDBInterface *instance();

  ~CDBInterface() override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  int UpdateRunNode(PHCompositeNode *topNode) override;

  void Disable() { disable = true; }
  void Enable() { disable = false; }

  void Disable_default() { disable_default = true; }
  void Enable_default() { disable_default = false; }

  std::string getUrl(const std::string &domain, const std::string &filename = "");

 private:
  CDBInterface(const std::string &name = "CDBInterface");

  static CDBInterface *__instance;
  SphenixClient *cdbclient{nullptr};
  bool disable{false};
  bool disable_default{false};
  std::set<std::tuple<std::string, std::string, uint64_t>> m_UrlVector;
};

#endif  // FFAMODULES_CDBINTERFACE_H
