#ifndef PHG4ParticleGeneratorBase_H__
#define PHG4ParticleGeneratorBase_H__

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <vector>

class PHG4InEvent;
class PHG4Particle;

class PHG4ParticleGeneratorBase: public SubsysReco
{
 public:
  virtual ~PHG4ParticleGeneratorBase();

  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);

  virtual void set_name(const std::string &particle = "proton");
  virtual void set_pid(const int pid);
  virtual void set_mom(const double x, const double y, const double z);
  virtual void set_vtx(const double x, const double y, const double z);
  virtual void set_t0(const double t) {t0 = t;}

  virtual void Print(const std::string &what = "ALL") const;
  virtual void AddParticle(const std::string &particle, const double x, const double y, const double z);
  virtual void AddParticle(const int pid, const double x, const double y, const double z);
  virtual void Embed(const int i=1) {embedflag = i;}
  void set_seed(const unsigned int iseed);
  unsigned int get_seed() const {return seed;}

 protected:
  PHG4ParticleGeneratorBase(const std::string &name="GENERATORBASE");
  int get_pdgcode(const std::string &name);
  std::string get_pdgname(const int pdgcode);
  void CheckAndCreateParticleVector();
  void SetParticleId(PHG4Particle *particle, PHG4InEvent *ineve);
  int embedflag;
  double vtx_x;
  double vtx_y;
  double vtx_z;
  double t0;
  std::vector<PHG4Particle *> particlelist;
  unsigned int seed;
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif
