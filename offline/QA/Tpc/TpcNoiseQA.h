// include guard to prevent TpcNoiseQA.h from being preprocessed multiple times
#ifndef QA_TPC_TPCNOISEQA_H
#define QA_TPC_TPCNOISEQA_H
//

// includes
#include <fun4all/SubsysReco.h>

#include <tpc/TpcMap.h>

#include <cmath>
#include <string>
#include <vector>

// Call classes to be used in code
class PHCompositeNode;
class TH2;
//

class TpcNoiseQA : public SubsysReco  // Inherit public parts of SubsysReco
{
  // list of public methods
 public:
  // Function that sets file name and allocated memory to adcSamples
  explicit TpcNoiseQA(const std::string &name = "TpcNoiseQA");

  ~TpcNoiseQA() override = default;

  // Initialization of run: where you fetch data and place histograms that
  // care about run number
  int InitRun(PHCompositeNode *topNode) override;

  // Where you do the desired work associated with the run (analysis)
  int process_event(PHCompositeNode *topNode) override;

  // called at the end of the run when processing is over
  int End(PHCompositeNode *topNode) override;

  // List of private members
 private:
  TpcMap M;

  TH2 *h_NPol_Ped_Mean{nullptr};  // Histogram of pedestal values on North Side
  TH2 *h_NPol_Ped_RMS{nullptr};   // Histogram of RMS values on North Side
  TH2 *h_SPol_Ped_Mean{nullptr};  // Histogram of pedestal values on South Side
  TH2 *h_SPol_Ped_RMS{nullptr};   // Histogram of RMS values on South Side

  const int r_bins_N{67};
// Number of radius bins
  double r_bins[68]{217.83, 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83, 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49, 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76, 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11, 1000};
  // Radius bin array
  const int nphi{205};
  // Number of phi bins
  double phi_bins[206]{0., 6.3083 - 2 * M_PI, 6.3401 - 2 * M_PI, 6.372 - 2 * M_PI, 6.4039 - 2 * M_PI, 6.4358 - 2 * M_PI, 6.4676 - 2 * M_PI, 6.4995 - 2 * M_PI, 6.5314 - 2 * M_PI, 0.2618, 0.2937, 0.3256, 0.3574, 0.3893, 0.4212, 0.453, 0.4849, 0.5168, 0.5487, 0.5805, 0.6124, 0.6443, 0.6762, 0.7081, 0.7399, 0.7718, 0.7854, 0.8173, 0.8491, 0.881, 0.9129, 0.9448, 0.9767, 1.0085, 1.0404, 1.0723, 1.1041, 1.136, 1.1679, 1.1998, 1.2317, 1.2635, 1.2954, 1.309, 1.3409, 1.3727, 1.4046, 1.4365, 1.4684, 1.5002, 1.5321, 1.564, 1.5959, 1.6277, 1.6596, 1.6915, 1.7234, 1.7552, 1.7871, 1.819, 1.8326, 1.8645, 1.8963, 1.9282, 1.9601, 1.992, 2.0238, 2.0557, 2.0876, 2.1195, 2.1513, 2.1832, 2.2151, 2.247, 2.2788, 2.3107, 2.3426, 2.3562, 2.3881, 2.42, 2.4518, 2.4837, 2.5156, 2.5474, 2.5793, 2.6112, 2.6431, 2.6749, 2.7068, 2.7387, 2.7706, 2.8024, 2.8343, 2.8662, 2.8798, 2.9117, 2.9436, 2.9754, 3.0073, 3.0392, 3.0711, 3.1029, 3.1348, 3.1667, 3.1986, 3.2304, 3.2623, 3.2942, 3.326, 3.3579, 3.3898, 3.4034, 3.4353, 3.4671, 3.499, 3.5309, 3.5628, 3.5946, 3.6265, 3.6584, 3.6903, 3.7221, 3.754, 3.7859, 3.8178, 3.8496, 3.8815, 3.9134, 3.927, 3.9589, 3.9907, 4.0226, 4.0545, 4.0864, 4.1182, 4.1501, 4.182, 4.2139, 4.2457, 4.2776, 4.3095, 4.3414, 4.3732, 4.4051, 4.437, 4.4506, 4.4825, 4.5143, 4.5462, 4.5781, 4.61, 4.6418, 4.6737, 4.7056, 4.7375, 4.7693, 4.8012, 4.8331, 4.865, 4.8968, 4.9287, 4.9606, 4.9742, 5.0061, 5.0379, 5.0698, 5.1017, 5.1336, 5.1654, 5.1973, 5.2292, 5.2611, 5.2929, 5.3248, 5.3567, 5.3886, 5.4204, 5.4523, 5.4842, 5.4978, 5.5297, 5.5615, 5.5934, 5.6253, 5.6572, 5.689, 5.7209, 5.7528, 5.7847, 5.8165, 5.8484, 5.8803, 5.9122, 5.944, 5.9759, 6.0078, 6.0214, 6.0533, 6.0851, 6.117, 6.1489, 6.1808, 6.2127, 6.2445, 6.2764, 2 * M_PI};  // Phi bin array

  double r_bins_new[2 * (67 + 1) + 205 + 1]{};

  int mod_arr[26]{2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};  // Array of TPC modules

  int FEE_map[26]{4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};  // FEE mapping for TPC modules

  float ave_adc_fee_channel[26][256]{};     // average ADC value in channel
  float std_adc_fee_channel[26][256]{};     // average RMS in channel
  float counts_adc_fee_channel[26][256]{};  // number of times channel appears

  bool dead{false};

  float temp1{0.0};
  float temp2{0.0};

  int feeM{0};
  unsigned int key{0};
  double R{0.0};
  double phi{0.0};
  double pedMean{0.0};
  double pedStdi{0.0};

  int side{0};                // Face of the TPC (0==North && 1==South)
  int m_nWaveformInFrame{0};  // Number of waveforms in a frame
  int m_Channel{0};           // Channel number
  int m_nSamples{0};          // Number of samples in waveform
  int m_FEE{0};

  void createHistos();
  std::string getHistoPrefix() const;

  std::string m_fname;  // Name of file given to program

  std::vector<unsigned short> m_adcSamples;  // ADC values in waveform
};

#endif  // ends preprocessing if TpcNoiseQA.h has already been defined
