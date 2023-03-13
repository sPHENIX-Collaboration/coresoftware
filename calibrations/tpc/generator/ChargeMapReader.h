class TH3;
template <class T>
class MultiArray;

//since we are never in the position of adding very large numbers to very small, floats are sufficient precision here.
class ChargeMapReader
{
 public:
  ChargeMapReader();  //calls the below with default values.
  ChargeMapReader(int _n0, float _rmin, float _rmax, int _n1, float _phimin, float _phimax, int _n2, float _zmin, float _zmax);
  ~ChargeMapReader();

 private:
  MultiArray<float>* charge = nullptr;
  TH3* hSourceCharge = nullptr;
  TH3* hChargeDensity = nullptr;
  float inputAxisScale = 1;    //multiply the r and z dimensions of the input histogram by this, when filling our internal array.  So if the input histogram is in mm and we want to fill our array in cm, inputUnit=0.1;
  float inputChargeScale = 1;  //multiply the content the input histogram bins by this, when filling our internal array.
  int nBins[3] = {1, 1, 1};    //r,phi,z bins of the output fixed-width array
  float lowerBound[3] = {0, 0, 0};
  float upperBound[3] = {999, 999, 999};
  float binWidth[3] = {999, 999, 999};

  bool CanInterpolateAt(float r, float phi, float z);  //checks whether it is okay to interpolate at this position in the charge density hist

  void RegenerateCharge();   //internal function to revise the internal array whenever the bounds change etc.
  void RegenerateDensity();  //internal function to rebuild the charge density map when the input map changes.

 public:
  static bool CanInterpolateAt(float x, float y, float z, TH3* h);  //checks whether it is okay to interpolate at this position in the supplied hist (a convenient utility)

  void FillChargeHistogram(TH3* h);  //fill the supplied histogram with the charge in the array.

  void AddChargeInBin(int r, int phi, int z, float q);
  void AddChargeAtPosition(float r, float phi, float z, float q);

  float GetChargeInBin(int r, int phi, int z);
  float GetChargeAtPosition(float r, float phi, float z);
  TH3* GetDensityHistogram() { return hChargeDensity; }  //returns the charge density hist if we still have it.

  bool ReadSourceCharge(const char* filename, const char* histname, float axisScale = 1., float contentScale = 1.);
  bool ReadSourceCharge(TH3* sourceHist, float axisScale = 1., float contentScale = 1.);

  void SetChargeInBin(int r, int phi, int z, float q);
  void SetChargeAtPosition(float r, float phi, float z, float q);

  bool SetOutputParameters(int _nr, float _rmin, float _rmax, int _nphi, float _phimin, float _phimax, int _nz, float _zmin, float _zmax);
  bool SetOutputBounds(float _rmin, float _rmax, float _phimin, float _phimax, float _zmin, float _zmax);
  bool SetOutputBins(int _nr, int _nphi, int _nz);
};
