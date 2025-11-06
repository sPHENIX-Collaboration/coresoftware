#ifndef JETBASE_TOWERJETINPUT_H
#define JETBASE_TOWERJETINPUT_H

#include "Jet.h"
#include "JetInput.h"

#include <calobase/RawTowerDefs.h>
#include <globalvertex/GlobalVertex.h>

#include <iostream>  // for cout, ostream

#include <vector>
// forward declarations
class PHCompositeNode;
class GlobalVertex;
class TowerJetInput : public JetInput
{
 public:
  TowerJetInput(Jet::SRC input, const std::string &prefix = "TOWERINFO_CALIB");
  ~TowerJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return m_input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

  void reset_GlobalVertexType()
  {
    m_use_vertextype = false;
    m_vertex_type.clear();
  }
  
  void set_GlobalVertexType(GlobalVertex::VTXTYPE type) 
  {
    reset_GlobalVertexType();
    m_use_vertextype = true;
    m_vertex_type.push_back(type);
  }

  void add_GlobalVertexType(GlobalVertex::VTXTYPE type)
  {
    m_use_vertextype = true;
    m_vertex_type.push_back(type);
  }

  void set_GlobalVertexTypes(const std::vector<GlobalVertex::VTXTYPE>& types)
  {
    reset_GlobalVertexType();
    m_use_vertextype = true;
    for(unsigned int i=0; i<types.size(); ++i)
    {
      m_vertex_type.push_back(types.at(i));
    }
  }

  float get_timing_e_threshold() { return m_timing_e_threshold; }
  
  void set_timing_e_threshold(float new_threshold) { m_timing_e_threshold = new_threshold; }

 private:
  Jet::SRC m_input;
  RawTowerDefs::CalorimeterId geocaloid{RawTowerDefs::CalorimeterId::NONE};
  bool m_use_towerinfo {false};
  std::string m_towerNodePrefix;
  std::string towerName;
  bool m_use_vertextype {false};
  std::vector<GlobalVertex::VTXTYPE> m_vertex_type{GlobalVertex::UNDEFINED};
  float m_timing_e_threshold{0.1};
};

#endif
