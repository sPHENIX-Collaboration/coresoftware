#ifndef G4JET_TRUTHJETINPUT_H
#define G4JET_TRUTHJETINPUT_H

#include "JetInput.h"

#include "Jet.h"

#include <iostream>  // for cout, ostream
#include <vector>

class PHCompositeNode;

class TruthJetInput : public JetInput
{
 public:
  TruthJetInput(Jet::SRC input);
  ~TruthJetInput() override {}

  //! by default, TruthJetInput process all truth primary particle.
  //! However, it can be configured to read only one or more embedded stream via add_embedding_flag()
  //! It can be useful for reconstruct truth jet for embedded pythia jets only, etc.
  //! Call add_embedding_flag() multiple times to add multiple embed stream
  void add_embedding_flag(const int embed_stream_id)
  {
    m_EmbedID.push_back(embed_stream_id);
  }

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return m_Input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

  void set_eta_range(float eta_min, float eta_max)
  {
    m_EtaMin = eta_min;
    m_EtaMax = eta_max;
  }

 private:
  Jet::SRC m_Input = Jet::VOID;
  float m_EtaMin = -4.;
  float m_EtaMax = 4.;

  //! if empty: process all primary particles
  //! if non-empty: only process primary particles in the selected embed stream.
  std::vector<int> m_EmbedID;

  bool use_embed_stream() { return m_EmbedID.size() > 0; }
};

#endif
