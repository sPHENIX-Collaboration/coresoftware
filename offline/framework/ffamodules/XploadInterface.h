// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_XPLOADINTERFACE_H
#define FFAMODULES_XPLOADINTERFACE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>  // for uint64_t
#include <set>
#include <string>
#include <tuple>  // for tuple

class PHCompositeNode;

class XploadInterface : public SubsysReco
{
 public:
  static XploadInterface *instance();

  ~XploadInterface() override {}

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  std::string getUrl(const std::string &domain, const std::string &filename = "");

 private:
  XploadInterface(const std::string &name = "XploadInterface");

  static XploadInterface *__instance;

  std::set<std::tuple<std::string, std::string, uint64_t>> m_UrlVector;
};

#endif  // FFAMODULES_XPLOADINTERFACE_H
