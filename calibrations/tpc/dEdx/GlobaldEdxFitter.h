#ifndef GLOBALDEDXFITTER_H
#define GLOBALDEDXFITTER_H

#include "bethe_bloch.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TChain.h"
#include "TGraph.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

class GlobaldEdxFitter
{
  public:
  GlobaldEdxFitter(double xmin = 10., double xmax = 50.) 
  {
    min_norm = xmin;
    max_norm = xmax;
  };
  void processResidualData(const std::string& infile, 
                  size_t ntracks = 200000, 
                  size_t skip = 0);
  void addTrack(double trk_dEdx, double trk_p);
  size_t getNtracks()
  {
    return dEdx.size();
  }

  double get_fitquality(double norm, double ZS_loss = 0.);
  double get_fitquality_new(double A);
  TF1* create_TF1(const std::string& name);
  TF2* create_TF2(const std::string& name);
  TF3* create_TF3_new(const std::string& name);
  double get_minimum();
  double get_minimum_new();
  std::pair<double,double> get_minimum_ZS();
  void set_range(double xmin, double xmax, double ZSmin, double ZSmax) 
  {
    min_norm = xmin;
    max_norm = xmax;
    min_ZS = ZSmin;
    max_ZS = ZSmax;
  }
  void reset()
  {
    p.clear();
    dEdx.clear();
  }
  std::vector<double> get_betagamma(double A);
  TGraph* graph_vsbetagamma(double A);
  TGraph* graph_vsp();
  private:
  std::vector<double> p;
  std::vector<double> dEdx;

  double get_fitquality_functor(const double* x);

  double get_fitquality_wrapper(double* x, double* par);
  double get_fitquality_wrapper_ZS(double* x, double* par);
  double get_fitquality_wrapper_new(double* x, double* par);
  double min_norm = 10.;
  double max_norm = 50.;
  double min_ZS = 0.;
  double max_ZS = 200.;
  double min_B = 8.;
  double max_B = 12.;
};

#endif // GLOBALDEDXFITTER_H
