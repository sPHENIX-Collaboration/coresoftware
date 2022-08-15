/**********************************/
/*  Program to arrange, init,     */
/*   and calc an MVA response     */
/*   Cameron Dean, LANL, 06/15/20 */
/**********************************/

#ifndef KFPARTICLESPHENIX_KFPARTICLEMVA_H
#define KFPARTICLESPHENIX_KFPARTICLEMVA_H

#include <Rtypes.h>  // for root's Float_t

#include <string>
#include <tuple>
#include <vector>

class KFParticle;
class KFPVertex;

namespace TMVA
{
  class Reader;
}

class KFParticle_MVA
{
 public:
  KFParticle_MVA() {}

  virtual ~KFParticle_MVA() {}

  std::tuple<TMVA::Reader *, std::vector<Float_t>> initMVA();

  Float_t evaluateMVA(TMVA::Reader *reader, std::vector<Float_t> reader_floats, KFParticle particle, KFPVertex vertex);

 protected:
  unsigned int m_nPars = 1;
  std::vector<std::string> m_mva_variable_list;
  std::string m_mva_type;
  std::string m_mva_path;

 private:
  unsigned int nMVApars = m_nPars;  //sizeof(m_mva_variable_list)/sizeof(m_mva_variable_list[0]);
  std::string method = m_mva_path + " method";
};

#endif  //KFPARTICLESPHENIX_KFPARTICLEMVA_H
