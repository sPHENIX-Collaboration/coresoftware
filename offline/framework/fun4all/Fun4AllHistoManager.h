// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLHISTOMANAGER_H
#define FUN4ALL_FUN4ALLHISTOMANAGER_H

#include "Fun4AllBase.h"

#include <map>
#include <string>

class TNamed;

class Fun4AllHistoManager : public Fun4AllBase
{
 public:
  explicit Fun4AllHistoManager(const std::string &name);
  ~Fun4AllHistoManager() override;

  void Print(const std::string &what = "ALL") const override;

  //! Register histogram or TTree object
  //! For histograms, enforce error calculation and propagation
  bool registerHisto(const std::string &hname, TNamed *h1d, const int replace = 0);

  //! Register histogram or TTree object
  //! For histograms, enforce error calculation and propagation
  bool registerHisto(TNamed *h1d, const int replace = 0);

  template <typename T>
  T *makeHisto(T *t)
  {
    if (not registerHisto(t))
    {
      delete t;
      t = nullptr;
    }
    return t;
  }
  int isHistoRegistered(const std::string &name) const;
  TNamed *getHisto(const std::string &hname) const;
  TNamed *getHisto(const unsigned int ihisto) const;
  std::string getHistoName(const unsigned int ihisto) const;
  unsigned int nHistos() const { return Histo.size(); }
  void Reset();
  int dumpHistos(const std::string &filename = "", const std::string &openmode = "RECREATE");
  void setOutfileName(const std::string &filename) { outfilename = filename; }

 private:
  std::string outfilename;
  std::map<const std::string, TNamed *> Histo;
};

#endif /* __FUN4ALLHISTOMANAGER_H */
