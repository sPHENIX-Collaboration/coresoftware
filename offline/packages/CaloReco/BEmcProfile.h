#include <vector>  // for vector

class EmcModule;
class TH1F;

class BEmcProfile
{

public:
  BEmcProfile(const char* fname);
  virtual ~BEmcProfile();

  float GetProb(std::vector<EmcModule>* plist, int NX, float en, float theta, float phi);
  float GetTowerEnergy( int iy, int iz, std::vector<EmcModule>* plist, int nx );
  void PredictEnergy(int ip, float en, float theta, float phi, float ddz, float ddy, float& ep, float& err);
  //  float GetProbTest(std::vector<EmcModule>* plist, int NX, float en, float theta, float& test_rr, float& test_et, float& test_ep, float& test_err);

protected:
  bool bloaded;

  float thresh;
  int nen;
  int nth;

  float* energy_array;
  float* theta_array;

  TH1F* *hmean;
  TH1F* *hsigma;
};
