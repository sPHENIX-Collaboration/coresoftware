    /**********************************/
   /*  Program to arrange, init,     */
  /*   and calc an MVA response     */
 /*   Cameron Dean, LANL, 06/15/20 */
/**********************************/

#ifndef KFParticle_MVA_H
#define KFParticle_MVA_H

//ROOT stuff
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

class KFParticle;
class KFParticleBase;
class KFPVertex;
//class TMVA;

class KFParticle_MVA 
{

public:

  KFParticle_MVA();

  ~KFParticle_MVA();

  tuple<TMVA::Reader*, vector<Float_t>> initMVA();

  Float_t evaluateMVA( TMVA::Reader *reader, vector<Float_t> reader_floats, KFParticle particle, KFPVertex vertex );

protected:

  unsigned int m_nPars;
  string m_mva_variable_list[ 99 ];
  string m_mva_type;
  string m_mva_path;

private:

  unsigned int nMVApars = m_nPars;  //sizeof(m_mva_variable_list)/sizeof(m_mva_variable_list[0]);
  string method = m_mva_path + " method";

};

#endif //KFParticle_MVA_H
