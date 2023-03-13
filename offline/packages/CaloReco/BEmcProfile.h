#include <string>
#include <vector>  // for vector

class EmcModule;
class TH1F;

class BEmcProfile
{
 public:
  explicit BEmcProfile(const std::string& fname);

  // delete copy ctor and assignment operator (cppcheck)
  explicit BEmcProfile(const BEmcProfile&) = delete;
  BEmcProfile& operator=(const BEmcProfile&) = delete;

  virtual ~BEmcProfile();

  float GetProb(std::vector<EmcModule>* plist, int NX, float en, float theta, float phi);
  float GetTowerEnergy(int iy, int iz, std::vector<EmcModule>* plist, int nx);
  void PredictEnergy(int ip, float en, float theta, float phi, float ddz, float ddy, float& ep, float& err);
  float PredictEnergyR(float energy, float theta, float phi, float rr);
  //  float GetProbTest(std::vector<EmcModule>* plist, int NX, float en, float theta, float& test_rr, float& test_et, float& test_ep, float& test_err);
  bool IsLoaded() { return bloaded; }
  int Verbosity() const { return m_Verbosity; }
  void Verbosity(const int i) { m_Verbosity = i; }

 protected:
  bool bloaded;

  float thresh;
  int nen;
  int nth;

  float* energy_array;
  float* theta_array;

  TH1F** hmean;
  TH1F** hsigma;
  TH1F** hr4;

 private:
  int m_Verbosity;
};
