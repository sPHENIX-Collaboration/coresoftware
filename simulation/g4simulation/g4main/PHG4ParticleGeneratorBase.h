// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGENERATORBASE_H
#define G4MAIN_PHG4PARTICLEGENERATORBASE_H

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_rng.h>

#include <string>  // for string
#include <vector>

class PHCompositeNode;
class PHG4InEvent;
class PHG4Particle;

class PHG4ParticleGeneratorBase : public SubsysReco
{
 public:
  ~PHG4ParticleGeneratorBase() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  virtual void set_name(const std::string &particle = "proton");
  virtual void set_pid(const int pid);
  virtual void set_mom(const double x, const double y, const double z);
  virtual void set_vtx(const double x, const double y, const double z);
  virtual void set_vtx_z(const double z) { m_Vtx_z = z; }
  virtual void set_t0(const double t) { m_TZero = t; }

  virtual double get_vtx_x() const { return m_Vtx_x; }
  virtual double get_vtx_y() const { return m_Vtx_y; }
  virtual double get_vtx_z() const { return m_Vtx_z; }
  virtual double get_t0() const { return m_TZero; }

  void Print(const std::string &what = "ALL") const override { PrintParticles(what); }
  virtual void PrintParticles(const std::string &what = "ALL") const;
  virtual void AddParticle(const std::string &particle, const double x, const double y, const double z);
  virtual void AddParticle(const int pid, const double x, const double y, const double z);
  virtual void Embed(const int i = 1) { m_EmbedFlag = i; }
  virtual int ReuseExistingVertex(PHCompositeNode *topNode);

  int get_reuse_existing_vertex() const { return m_ReUseExistingVertexFlag; }
  void set_reuse_existing_vertex(const int i = 1) { m_ReUseExistingVertexFlag = i; }
  void set_seed(const unsigned int iseed);
  unsigned int get_seed() const { return m_Seed; }

 protected:
  PHG4ParticleGeneratorBase(const std::string &name = "GENERATORBASE");
  int get_pdgcode(const std::string &name) const;
  std::string get_pdgname(const int pdgcode) const;
  double get_mass(const int pdgcode) const;
  void CheckAndCreateParticleVector();
  void SetParticleId(PHG4Particle *particle, PHG4InEvent *ineve);
  gsl_rng *RandomGenerator() const { return m_RandomGenerator; }
  int EmbedFlag() const { return m_EmbedFlag; }
  std::vector<PHG4Particle *>::iterator particlelist_begin() { return particlelist.begin(); }
  std::vector<PHG4Particle *>::iterator particlelist_end() { return particlelist.end(); }
  void ResetParticleList();

 private:
  gsl_rng *m_RandomGenerator = nullptr;
  int m_EmbedFlag = 0;
  int m_ReUseExistingVertexFlag = 0;
  unsigned int m_Seed = 0;
  double m_Vtx_x = 0.;
  double m_Vtx_y = 0.;
  double m_Vtx_z = 0.;
  double m_TZero = 0.;
  std::vector<PHG4Particle *> particlelist;
};

#endif
