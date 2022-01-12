#ifndef G4EVAL_DSTCOMPRESSOR_H
#define G4EVAL_DSTCOMPRESSOR_H

#include "compressor_generator.h"
#include "RtypesCore.h"


class DSTCompressor {
public:
  // whether to print the debug info
  const bool debug = false;
  // exit prematurely after this much residuals
  const int earlybreak = 1000 * 1000 * 1000;

  const float dist_phi = 0.05;
  const float dist_z = 3;
  int phiBits;
  int zBits;



  using bit5 = std::bitset<5>;
  using bit8 = std::bitset<8>;
  using bit10 = std::bitset<10>;
  using bstream = std::vector<u_char>;
  using Dict = std::vector<float>;

  const int wordwidth = 8; // a char has 8 bits
  // phi mean, phi sigma, z mean, z sigma
  using Pars = std::array<double, 4>;

  Dict phiDict;
  Dict zDict;

  DSTCompressor(float phiMean = -3e-5,
                float phiSigma = 0.015,
                float zMean = -0.00063,
                float zSigma = 0.07036,
                int nBitPhi = 8,
                int nBitZ = 8)
  : phiBits(nBitPhi),
    zBits(nBitZ)
{
  int numPoints = 1000000;
  phiDict = compress_gaussian_dist(phiMean, phiSigma, numPoints, phiBits);
  zDict = compress_gaussian_dist(zMean, zSigma, numPoints, zBits);
}

/** Generate the lookup table with a Gaussian model
 */
Dict compress_gaussian_dist(double mean, double stddev,
                                            int numPoints, int numBits) {
  std::vector<ushort> order;
  std::vector<float> dict;
  std::vector<unsigned long> cnt;
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(mean, stddev);
  float maxAbsError = approx(&order, &dict, &cnt, numPoints, generator,
                               distribution, (size_t)pow(2, numBits));
  std::cout << "Compressing with " << numBits << " bits" << std::endl;
  std::cout << "Number of clusters = " << dict.size() << std::endl;
  std::cout << "Maximum absolute error = " << maxAbsError << std::endl;
  return dict;
}

  unsigned short compressPhi(float inPhi) {
    return residesIn(inPhi, &phiDict);
  };

  unsigned short compressZ(float inZ) {
    return residesIn(inZ, &zDict);
  };

  float decompressPhi(unsigned short key) {
    return phiDict.at(key);
  }

  float decompressZ(unsigned short key) {
    return zDict.at(key);
  }

};

#endif // G4EVAL_DSTCOMPRESSOR_H
