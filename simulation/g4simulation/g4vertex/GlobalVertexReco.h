// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEXRECO_H
#define G4VERTEX_GLOBALVERTEXRECO_H

//===========================================================
/// \file GlobalVertexReco.h
/// \brief reconstruct the best possible vertexes by
/// combining results from multiple detectors
/// \author Mike McCumber
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

/// \class GlobalVertexReco
///
/// \brief simple truth vertex smearing algorithm
///
class GlobalVertexReco : public SubsysReco
{
 public:
  GlobalVertexReco(const std::string &name = "GlobalVertexReco");
  ~GlobalVertexReco() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

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

  float _xdefault, _xerr;
  float _ydefault, _yerr;
  float _tdefault, _terr;
};

#endif  // G4VERTEX_GLOBALVERTEXRECO_H
