#ifndef FLOWAFTERBURNER_FLOWAFTERBURNER_H
#define FLOWAFTERBURNER_FLOWAFTERBURNER_H

#include <string>
#include <array>

#include "AfterburnerAlgo.h"

namespace CLHEP
{
  class HepRandomEngine;
}
namespace HepMC
{
  class GenEvent;
  class GenParticle;
}

class Afterburner
{

 public:

  explicit  Afterburner( const std::string &algorithmName = "MINBIAS",
              CLHEP::HepRandomEngine *engine = nullptr,
              float mineta = -5.0f, float maxeta = 5.0f,
              float minpt = 0.0f, float maxpt = 100.0f );

  ~Afterburner();
  
  Afterburner(const Afterburner&) = delete;
  Afterburner& operator=(const Afterburner&) = delete;
  Afterburner(Afterburner&&) noexcept;
  Afterburner& operator=(Afterburner&&) noexcept;


  // overloaded algo setter (legacy)
  void setAlgo(AfterburnerAlgo::flowAfterburnerAlgorithm algo_type);
  void setAlgo(const std::string &name);
  void setAlgo(AfterburnerAlgo * algo);
  
  AfterburnerAlgo * getAlgo() { return m_algo; }
  
  void setEngine(CLHEP::HepRandomEngine *engine);
  CLHEP::HepRandomEngine * getEngine() { return m_engine; }

  void setEtaRange(float mineta, float maxeta);
  void setPtRange(float minpt, float maxpt);

  float getPsiN(unsigned int n) const;
 
  static double vn_func(double x, void *params);
  void throw_psi_n(HepMC::GenEvent *event);

  void AddFlowToParentAndMoveDescendants(HepMC::GenEvent *event, HepMC::GenParticle *parent);

  // add all the extra arguments for legacy flowAfterburner function
  int flowAfterburner(HepMC::GenEvent *event, 
                      CLHEP::HepRandomEngine *engine = nullptr,
                      const std::string &algorithmName = "",
                      float mineta = -5.0f, float maxeta = 5.0f,
                      float minpt = 0.0f, float maxpt = 100.0f);

 private:
    
  AfterburnerAlgo * m_algo = nullptr;
  CLHEP::HepRandomEngine * m_engine = nullptr;
  bool m_ownAlgo = false;
  bool m_ownEngine = false;
  float m_mineta = -5.0f;
  float m_maxeta = 5.0f;
  float m_minpt = 0.0f;
  float m_maxpt = 100.0f;
  double m_phishift = 0.0; // shift of the reaction plane angle in phi, used to align with the impact parameter

  void setPsiN(unsigned int n, float psi);
  float m_psi_n[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; // reaction plane angles

  // Legacy arguments
  void readLegacyArguments(
      CLHEP::HepRandomEngine *engine,
      const std::string &algorithmName,
      float mineta, float maxeta,
      float minpt, float maxpt);  
};

// legacy function, use the Afterburner class instead
int flowAfterburner(HepMC::GenEvent *inEvent,
                    CLHEP::HepRandomEngine *engine,
                    const std::string &algorithmName,
                    float mineta, float maxeta,
                    float minpt, float maxpt);

#endif
