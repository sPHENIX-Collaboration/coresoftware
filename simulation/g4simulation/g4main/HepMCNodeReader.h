#ifndef HEPMCNODEREADER_H__
#define HEPMCNODEREADER_H__

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <string>

class PHCompositeNode;

class HepMCNodeReader : public SubsysReco
{
 public:
  HepMCNodeReader(const std::string &name = "HEPMCREADER");
  virtual ~HepMCNodeReader();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  void Embed(const int i=1) {_embed_flag = i;}
  void VertexPosition(const double v_x, const double v_y, const double v_z);


  void SmearVertex(const double s_x, const double s_y, const double s_z);

private:
  double smeargauss(const double width);
  double smearflat(const double width);
  int _embed_flag;
  double vertex_pos_x;
  double vertex_pos_y;
  double vertex_pos_z;
  double width_vx;
  double width_vy;
  double width_vz;

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
 
};

#endif

