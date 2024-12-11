// header file for biquad filter. Focus will be on low pass filter but code will allow for other types of filters if needed

#ifndef Biquad_h
#define Biquad_h

// enumerated types for each filter

enum {
    bq_type_lowpass = 0,
    bq_type_highpass,
    bq_type_bandpass,
    bq_type_notch,
    bq_type_peak,
    bq_type_lowshelf,
    bq_type_highshelf
};

// Definition of Biquad class and member functions

class Biquad {
 public:
  Biquad(); // Bare-bones constructor of Biquad class
  Biquad(int ttype, double tFc, double tQ, double tpeakGainDB); // Class constructor with tunable input parameters. Useful if parameters are fixed
  ~Biquad(); // Destructor of Biquad class
  void setType(int ttype); // Function to set filter type
  void setQ(double tQ); // Function to set frequency range of filter. Smaller Q give wider range, larger Q gives narrower
  void setFc(double tFc); // Function to set normalized cutoff frequency (Fc = cutoff/sample rate)
  void setPeakGain(double tpeakGainDB); // Function to set gain at filtering point of filter
  void setBiquad(int ttype, double tFc, double tQ, double tpeakGainDB); // Function to set parameters of filter. Useful if parameters aren't fixed
  float process(float in) ; // Function to run input values through filter
  
 protected:
  void calcBiquad(void); // Function to set initial variables of each filter type
  
  // Initial variables for filters
  int type;
  double a0, a1, a2, b1, b2;
  double Fc, Q, peakGain;
  double z1, z2;
};

// Inline function for processing data through the filter. Should be small/fast running function so inline prevents overhead

inline float Biquad::process(float in) {
  double out = in * a0 + z1; // Filter output
  z1 = in * a1 + z2 - b1 * out; // updated parameter z1 using output
  z2 = in * a2 - b2 * out; // updated parameter z2 using output
  return out; // return output of the filtering
}

#endif // Biquad_h
