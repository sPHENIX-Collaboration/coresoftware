/**********************************/
/*  Program to arrange, init,     */
/*   and calc an MVA response     */
/*   Cameron Dean, LANL, 06/15/20 */
/**********************************/
#include "KFParticle_MVA.h"
#include "KFParticle_Tools.h"

#include <KFPVertex.h>
#include <KFParticle.h>

#include <TMVA/Reader.h>  // for Reader
#include <TMVA/Tools.h>   // for Tools

#include <algorithm>  // for max
#include <map>
#include <memory>   // for allocator_traits<>::value_type
#include <utility>  // for pair

KFParticle_Tools kfpTools;

std::tuple<TMVA::Reader *, std::vector<Float_t>> KFParticle_MVA::initMVA()
{
  TMVA::Tools::Instance();  //Start TMVA
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

  std::vector<Float_t> reader_floats;

  for (unsigned int i = 0; i < nMVApars; ++i)
  {
    reader_floats.push_back(0);
    reader->AddVariable(m_mva_variable_list[i].c_str(), &reader_floats[i]);
  }
  reader->BookMVA(method.c_str(), m_mva_path.c_str());

  return make_tuple(reader, reader_floats);
}

Float_t KFParticle_MVA::evaluateMVA(TMVA::Reader *reader, std::vector<Float_t> reader_floats, KFParticle particle, KFPVertex vertex)
{
  KFParticle kfpvertex(vertex);
  std::map<std::string, float> possibleVariables =
      {
          {"motherIPchi2", particle.GetDeviationFromVertex(kfpvertex)},
          {"motherFDchi2", kfpTools.flightDistanceChi2(particle, vertex)}};

  for (unsigned int iPar = 0; iPar < nMVApars; ++iPar) reader_floats[iPar] = possibleVariables.find(m_mva_variable_list[iPar])->second;

  return (Float_t) reader->EvaluateMVA(method.c_str());
}
