#include <fastjet/ClusterSequence.hh>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace fastjet;

bool enablep = false;
struct algo_info
{
  JetAlgorithm algorithm;
  double R;
  int PID;
  double EtaMin;
  double EtaMax;
  double EtMin;
};
vector<algo_info> algo_info_vec;

std::map<std::string, JetAlgorithm> algorithms;

struct loaderObj
{
  loaderObj() {
    static bool init = false;
    if (!init)
      {
	algorithms["KT"] = kt_algorithm;
	algorithms["CAMBRIDGE"] = cambridge_algorithm;
	algorithms["ANTIKT"] = antikt_algorithm;
	algorithms["GENKT"] = genkt_algorithm;
	algorithms["CAMBRIDGE_FOR_PASSIVE"] = cambridge_for_passive_algorithm;
	algorithms["GENKT_FOR_PASSIVE"] = genkt_for_passive_algorithm;
	algorithms["EE_KT"] = ee_kt_algorithm;
	algorithms["EE_GENKT"] = ee_genkt_algorithm;
	algorithms["PLUGIN"] = plugin_algorithm ;
      }
  }
};

loaderObj loader;

void
hijfst_control(int enable, vector<string> valgorithm, vector<float> vR, vector<int> vPID, vector<float> vEtaMin, vector<float> vEtaMax, vector<float> vEtMin)
{
  enablep = (enable==1) ? true: false;
  
  algo_info_vec.clear();
  for (unsigned int i = 0; i < valgorithm.size(); ++i)
    {
      string algorithmName = valgorithm[i];
      transform(algorithmName.begin(), algorithmName.end(), algorithmName.begin(), ::toupper);
      algo_info algo;
      algo.algorithm = ((algorithms.find(algorithmName) == algorithms.end()) ? antikt_algorithm : algorithms[algorithmName]);
      algo.R = vR[i];
      algo.PID = vPID[i];
      algo.EtaMin = vEtaMin[i];
      algo.EtaMax = vEtaMax[i];
      algo.EtMin = vEtMin[i];
      algo_info_vec.push_back(algo);
    }
}

extern "C" 
void
hijfst_(int *n, int *N, int *K, float *P, float *V)
{
  if (!enablep) return;
  const int M = *N;
    
  int *K1 = new(K)     int[M];
  int *K2 = new(K+M)   int[M];
  int *K3 = new(K+2*M) int[M];
  int *K4 = new(K+3*M) int[M];
  int *K5 = new(K+4*M) int[M];

  float *px = new(P)     float[M];
  float *py = new(P+M)   float[M];
  float *pz = new(P+2*M) float[M];
  float *E  = new(P+3*M) float[M];
  float *m  = new(P+4*M) float[M];

  float *V1 = new(V)     float[M];
  float *V2 = new(V+M)   float[M];
  float *V3 = new(V+2*M) float[M];
  float *V4 = new(V+3*M) float[M];
  float *V5 = new(V+4*M) float[M];

  std::vector<PseudoJet> particles;
    
  for (int i = 0; i < *n; i++)
    {
      // We only want real, final state particles
      if (K1[i] == 1)
	{
	  particles.push_back(PseudoJet(px[i], py[i], pz[i], E[i]));
	}
    }

  for (vector<algo_info>::iterator it = algo_info_vec.begin(); it != algo_info_vec.end(); ++it)
    {
      algo_info a = *it;

      // Set the algorithm and R value
      JetDefinition jet_def(a.algorithm, a.R);
      
      ClusterSequence cs(particles, jet_def);
      std::vector<PseudoJet> jets = cs.inclusive_jets();
    
      // Loop over all the "jets" and add their kinematic properties to
      // the LUJET common block.  Set their status to 103 to distinguish
      // them.
      for (unsigned i = 0; i < jets.size(); i++) 
	{
	  // These cuts should be configurable!
//        if (abs(jets[i].eta()) > 2.0 or jets[i].E() < 5.0)
      if (jets[i].eta() < a.EtaMin or jets[i].eta() > a.EtaMax or  jets[i].Et() < a.EtMin)
	    {
	      continue;
	    }
	  int j = (*n)++;
	  
	  K1[j] = 103;
	  K2[j] = a.PID;
	  K3[j] = 0;
	  K4[j] = 0;
	  K5[j] = 0;

	  px[j] = jets[i].px();
	  py[j] = jets[i].py();
	  pz[j] = jets[i].pz();
	  E[j]  = jets[i].E();
	  m[j]  = jets[i].m();

	  V1[j] = 0.0;
	  V2[j] = 0.0;
	  V3[j] = 0.0;
	  V4[j] = 0.0;
	  V5[j] = 0.0;
	}    
    }
}
  
