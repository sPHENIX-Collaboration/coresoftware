// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_GLOBALVERTEXRECO_H
#define GLOBALVERTEX_GLOBALVERTEXRECO_H

//===========================================================
/// \file GlobalVertexReco.h
/// \brief reconstruct the best possible vertexes by
/// combining results from multiple detectors
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>  // for string

class PHCompositeNode;

/// \class GlobalVertexReco
///
/// \brief simple truth vertex smearing algorithm
///
class GlobalVertexReco : public SubsysReco
{
 public:
  GlobalVertexReco(const std::string &name = "GlobalVertexReco");
  ~GlobalVertexReco() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_x_defaults(float xdefault, float xerr)
  {
    _xdefault = xdefault;
    _xerr = xerr;
  }
  void set_y_defaults(float ydefault, float yerr)
  {
    _ydefault = ydefault;
    _yerr = yerr;
  }
  void set_t_defaults(float tdefault, float terr)
  {
    _tdefault = tdefault;
    _terr = terr;
  }

 private:
  int CreateNodes(PHCompositeNode *topNode);

  float _xdefault{0.};
  float _xerr{0.3};
  float _ydefault{0.};
  float _yerr{0.3};
  float _tdefault{0.};
  float _terr{0.2};
};

#endif  // GLOBALVERTEX_GLOBALVERTEXRECO_H
