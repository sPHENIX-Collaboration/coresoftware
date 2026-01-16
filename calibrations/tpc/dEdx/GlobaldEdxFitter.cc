#include "GlobaldEdxFitter.h"

#include "bethe_bloch.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TChain.h"
#include "TGraph.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

void GlobaldEdxFitter::processResidualData(const std::string& infile, size_t ntracks, size_t skip)
{
  std::unique_ptr<TChain> t = std::make_unique<TChain>();
  t->Add((infile+"?#residualtree").c_str());
//  TFile* f = TFile::Open(infile.c_str());
//  TTree* t = (TTree*)f->Get("residualtree");

  float px;
  float py;
  float pz;
  float dedx;
  float eta;
  int nmaps;
  int nintt;
  int ntpc;
  float dcaxy;

  t->SetBranchAddress("px",&px);
  t->SetBranchAddress("py",&py);
  t->SetBranchAddress("pz",&pz);
  t->SetBranchAddress("dedx",&dedx);
  t->SetBranchAddress("eta",&eta);
  t->SetBranchAddress("nmaps",&nmaps);
  t->SetBranchAddress("nintt",&nintt);
  t->SetBranchAddress("ntpc",&ntpc);
  t->SetBranchAddress("dcaxy",&dcaxy);

  size_t total_entries = t->GetEntriesFast();

  for(size_t entry=skip; entry<(skip+ntracks); entry++)
  {
    if(entry==total_entries)
    {
      break;
    }
    if(entry % 1000 == 0)
    {
      std::cout << entry << std::endl;
    }
    t->GetEntry(entry);
    if(nmaps>0 && nintt>0 && fabs(eta)<1. && dcaxy<0.5 && ntpc>30)
    {
      p.push_back(sqrt(px*px+py*py+pz*pz));
      dEdx.push_back(dedx);
    }
  }
  std::cout << "number of good tracks: " << p.size() << std::endl;
  //f->Close();
}

void GlobaldEdxFitter::addTrack(double trk_dEdx, double trk_p)
{
  dEdx.push_back(trk_dEdx);
  p.push_back(trk_p);
}

double GlobaldEdxFitter::get_fitquality_new(double A)
{
  //double chi2 = 0.;
  //double ndf = -1.;

  double pi_chi2 = 0.;
  double K_chi2 = 0.;
  double p_chi2 = 0.;
  double d_chi2 = 0.;
  double pi_ndf = -1.;
  double K_ndf = -1.;
  double p_ndf = -1.;
  double d_ndf = -1.;

  for(size_t i=0; i<dEdx.size(); i++)
  {
    const double dedx_pi = bethe_bloch_new_1D(p[i]/dedx_constants::m_pi,A);
    const double dedx_K = bethe_bloch_new_1D(p[i]/dedx_constants::m_K,A);
    const double dedx_p = bethe_bloch_new_1D(p[i]/dedx_constants::m_p,A);
    const double dedx_d = bethe_bloch_new_1D(p[i]/dedx_constants::m_d,A);

    //std::cout << "dedx: (" << dedx_pi << ", " << dedx_K << ", " << dedx_p << ", " << dedx_d << ")" << std::endl;
    //std::cout << "measured: " << dEdx[i] << std::endl;

    const double pi_dist = fabs(dEdx[i]-dedx_pi);
    const double K_dist = fabs(dEdx[i]-dedx_K);
    const double p_dist = fabs(dEdx[i]-dedx_p);
    const double d_dist = fabs(dEdx[i]-dedx_d);

    //std::cout << "dist: (" << pi_dist << ", " << K_dist << ", " << p_dist << ", " << d_dist << ")" << std::endl;

    if(pi_dist<K_dist && pi_dist<p_dist && pi_dist<d_dist)
    {
      //std::cout << "pion" << std::endl;
      //chi2 += pi_dist*pi_dist;
      pi_chi2 += pi_dist*pi_dist;
      pi_ndf += 1.;
    }
    else if(K_dist<pi_dist && K_dist<p_dist && K_dist<d_dist)
    {
      //std::cout << "K" << std::endl;
      //chi2 += K_dist*K_dist;
      K_chi2 += K_dist*K_dist;
      K_ndf += 1.;
    }
    else if(p_dist<pi_dist && p_dist<K_dist && p_dist<d_dist)
    {
      //std::cout << "p" << std::endl;
      //chi2 += p_dist*p_dist;
      p_chi2 += p_dist*p_dist;
      p_ndf += 1.;
    }

    else if(d_dist<pi_dist && d_dist<K_dist && d_dist<p_dist)
    {
      //std::cout << "d" << std::endl;
      //chi2 += d_dist*d_dist;
      d_chi2 += d_dist*d_dist;
      d_ndf += 1.;
    }

    //ndf += 1.;
    //std::cout << "chi2: " << chi2 << std::endl;
    //std::cout << "ndf: " << ndf << std::endl;
  }

  if(pi_ndf<1.) {pi_chi2 = 0.; pi_ndf = 1.;}
  if(K_ndf<1.) {K_chi2 = 0.; K_ndf = 1.;}
  if(p_ndf<1.) {p_chi2 = 0.; p_ndf = 1.;}
  if(d_ndf<1.) {d_chi2 = 0.; d_ndf = 1.;}

  const double quality = sqrt((pi_chi2*pi_chi2)/(pi_ndf*pi_ndf)+(K_chi2*K_chi2)/(K_ndf*K_ndf)+(p_chi2*p_chi2)/(p_ndf*p_ndf)+(d_chi2*d_chi2)/(d_ndf*d_ndf));
  //const double quality = chi2/ndf;

  //std::cout << "A: " << A << " quality: " << quality << std::endl;

  return quality;
}

double GlobaldEdxFitter::get_fitquality(double norm, double ZS_loss)
{
  float pi_chi2 = 0.;
  float K_chi2 = 0.;
  float p_chi2 = 0.;
  float d_chi2 = 0.;
  float pi_ndf = -1.;
  float K_ndf = -1.;
  float p_ndf = -1.;
  float d_ndf = -1.;
  //float ndf = -1;

  for(size_t i=0; i<dEdx.size(); i++)
  {

    const float dedx_pi  = norm*bethe_bloch_total(p[i]/dedx_constants::m_pi) - ZS_loss;
    const float dedx_K = norm*bethe_bloch_total(p[i]/dedx_constants::m_K) - ZS_loss;
    const float dedx_p = norm*bethe_bloch_total(p[i]/dedx_constants::m_p) - ZS_loss;
    const float dedx_d = norm*bethe_bloch_total(p[i]/dedx_constants::m_d) - ZS_loss;

    const float pi_dist = fabs(dedx_pi-dEdx[i]);
    const float K_dist = fabs(dedx_K-dEdx[i]);
    const float p_dist = fabs(dedx_p-dEdx[i]);
    const float d_dist = fabs(dedx_d-dEdx[i]);

    if(pi_dist<K_dist && pi_dist<p_dist && pi_dist<d_dist)
    {
      pi_chi2 += pi_dist*pi_dist/fabs(dedx_pi);
      pi_ndf += 2.;
    }
    else if(K_dist<pi_dist && K_dist<p_dist && K_dist<d_dist)
    {
      K_chi2 += K_dist*K_dist/fabs(dedx_K);
      K_ndf += 2.;
    }
    else if(p_dist<pi_dist && p_dist<K_dist && p_dist<d_dist)
    {
      p_chi2 += p_dist*p_dist/fabs(dedx_p);
      p_ndf += 2.;
    }
    else if(d_dist<pi_dist && d_dist<K_dist && d_dist<p_dist)
    {
      d_chi2 += d_dist*d_dist/fabs(dedx_d);
      d_ndf += 2.;
    }
    //ndf += 2;
  }

  const double quality = sqrt((pi_chi2*pi_chi2)/(pi_ndf*pi_ndf) +
                        (K_chi2*K_chi2)/(K_ndf*K_ndf) +
                        (p_chi2*p_chi2)/(p_ndf*p_ndf) +
                        (d_chi2*d_chi2)/(d_ndf*d_ndf));

/*
  const double quality = sqrt((pi_chi2*pi_chi2 +
                        K_chi2*K_chi2 +
                        p_chi2*p_chi2)
                        /(ndf*ndf));
*/
  //std::cout << "norm: " << norm << " ZS: " << ZS_loss << std::endl;
  //std::cout << "quality: " << quality << std::endl;

  return quality;
}

std::vector<double> GlobaldEdxFitter::get_betagamma(double A)
{
  std::vector<double> betagamma;
  for(size_t i=0; i<dEdx.size(); i++)
  {
    const double dedx_pi = bethe_bloch_new_1D(p[i]/dedx_constants::m_pi,A);
    const double dedx_K = bethe_bloch_new_1D(p[i]/dedx_constants::m_K,A);
    const double dedx_p = bethe_bloch_new_1D(p[i]/dedx_constants::m_p,A);
    const double dedx_d = bethe_bloch_new_1D(p[i]/dedx_constants::m_d,A);

    //std::cout << "dedx: (" << dedx_pi << ", " << dedx_K << ", " << dedx_p << ", " << dedx_d << ")" << std::endl;
    //std::cout << "measured: " << dEdx[i] << std::endl;

    const double pi_dist = fabs(dEdx[i]-dedx_pi);
    const double K_dist = fabs(dEdx[i]-dedx_K);
    const double p_dist = fabs(dEdx[i]-dedx_p);
    const double d_dist = fabs(dEdx[i]-dedx_d);

    if(pi_dist<K_dist && pi_dist<p_dist && pi_dist<d_dist)
    {
      betagamma.push_back(p[i]/dedx_constants::m_pi);
    }
    else if(K_dist<pi_dist && K_dist<p_dist && K_dist<d_dist)
    {
      betagamma.push_back(p[i]/dedx_constants::m_K);
    }
    else if(p_dist<pi_dist && p_dist<K_dist && p_dist<d_dist)
    {
      betagamma.push_back(p[i]/dedx_constants::m_p);
    }
    else if(d_dist<pi_dist && d_dist<K_dist && d_dist<p_dist)
    {
      betagamma.push_back(p[i]/dedx_constants::m_d);
    }
  }
  return betagamma;
}

double GlobaldEdxFitter::get_fitquality_functor(const double* x)
{
  return get_fitquality_new(x[0]);
}

double GlobaldEdxFitter::get_fitquality_wrapper(double* x, double* par)
{
  return get_fitquality(x[0]);
}

double GlobaldEdxFitter::get_fitquality_wrapper_ZS(double* x, double* par)
{
  return get_fitquality(x[0],x[1]);
}

double GlobaldEdxFitter::get_fitquality_wrapper_new(double* x, double* par)
{
  return get_fitquality_new(x[0]);
}

TF3* GlobaldEdxFitter::create_TF3_new(const std::string& name)
{
  TF3* f = new TF3(name.c_str(),this,&GlobaldEdxFitter::get_fitquality_wrapper_new,min_norm,max_norm,min_B,max_B,min_ZS,max_ZS,0,"GlobaldEdxFitter","get_fitquality_wrapper_new");
  return f;
}

TF1* GlobaldEdxFitter::create_TF1(const std::string& name)
{
  TF1* f = new TF1(name.c_str(),this,&GlobaldEdxFitter::get_fitquality_wrapper_new,min_norm,max_norm,0,"GlobaldEdxFitter","get_fitquality_wrapper");
  return f;
}

TF2* GlobaldEdxFitter::create_TF2(const std::string& name)
{
  TF2* f = new TF2(name.c_str(),this,&GlobaldEdxFitter::get_fitquality_wrapper_ZS,min_norm,max_norm,min_ZS,max_ZS,0,"GlobaldEdxFitter","get_fitquality_wrapper_ZS");
  return f;
}

double GlobaldEdxFitter::get_minimum_new()
{
/*
  TF3* f = create_TF3_new("temp");
  double minA;
  double minB;
  double minC;
  f->GetMinimumXYZ(minA,minB,minC);
  delete f;
  return std::make_tuple(minA,minB,minC);
*/
  ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.1);
  minimizer->SetPrintLevel(1);
  ROOT::Math::Functor f(this,&GlobaldEdxFitter::get_fitquality_functor,1);
  double step[1] = {.01};
  double variable[1] = {20.};
  minimizer->SetFunction(f);
  minimizer->SetVariable(0,"A",variable[0],step[0]);
  minimizer->Minimize();
  const double *xs = minimizer->X();
  delete minimizer;
  return xs[0];
}

double GlobaldEdxFitter::get_minimum()
{
  TF1* f = create_TF1("temp");
  f->SetNpx(1000);
  double minX = f->GetMinimumX();
  delete f;
  return minX;
}

std::pair<double,double> GlobaldEdxFitter::get_minimum_ZS()
{
  TF2* f = create_TF2("temp");
  double minX;
  double minY;
  f->GetMinimumXY(minX,minY);
  delete f;
  return std::make_pair(minX,minY);
}

TGraph* GlobaldEdxFitter::graph_vsbetagamma(double A)
{
  std::vector<double> betagamma = get_betagamma(A);
  TGraph* g = new TGraph(dEdx.size(),betagamma.data(),dEdx.data());
  return g;
}

TGraph* GlobaldEdxFitter::graph_vsp()
{
  TGraph* g = new TGraph(dEdx.size(),p.data(),dEdx.data());
  return g;
}
