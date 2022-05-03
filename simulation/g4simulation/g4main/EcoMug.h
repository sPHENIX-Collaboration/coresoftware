/////////////////////////////////////////////////////////////////////////////////////
// EcoMug: Efficient COsmic MUon Generator                                         //
// Copyright (C) 2021 Davide Pagano <davide.pagano@unibs.it>                       //
// EcoMug is based on the following work:                                          //
// D. Pagano, G. Bonomi, A. Donzella, A. Zenoni, G. Zumerle, N. Zurlo,             //
// "EcoMug: an Efficient COsmic MUon Generator for cosmic-ray muons applications", //
// doi:10.1016/j.nima.2021.165732                                                  //
//                                                                                 //
// This program is free software: you can redistribute it and/or modify            //
// it under the terms of the GNU General Public License as published by            //
// the Free Software Foundation, either version 3 of the License, or               //
// (at your option) any later version.                                             //
//                                                                                 //
// This program is distributed in the hope that it will be useful,                 //
// but WITHOUT ANY WARRANTY; without even the implied warranty of                  //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   //
// GNU General Public License for more details.                                    //
//                                                                                 //
// You should have received a copy of the GNU General Public License               //
// along with this program.  If not, see <https://www.gnu.org/licenses/>.          //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef EcoMug_H
#define EcoMug_H

//#include <math.h>
#include <array>
#include <functional>
#include <random>

//! Fast generation of random numbers
//! This class is based on the xoroshiro128+ generator.
//! https://prng.di.unimi.it/
class EMRandom
{
 public:
  EMRandom()
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, std::numeric_limits<uint64_t>::max());
    s[0] = 12345.;  //dis(gen);
    s[1] = 12345;   //dis(gen);
  };

  void SetSeed(uint64_t seed)
  {
    s[0] = seed;
    s[1] = seed;
  };

  double GenerateRandomDouble()
  {
    uint64_t x = next();
    return to_double(x);
  };

  double GenerateRandomDouble(double x1, double x2)
  {
    return (x2 - x1) * GenerateRandomDouble() + x1;
  };

  int64_t rotl(const uint64_t x, int k)
  {
    return (x << k) | (x >> (64 - k));
  };

  uint64_t next()
  {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = rotl(s1, 36);
    return result;
  };

  double to_double(uint64_t x)
  {
    union U
    {
      uint64_t i;
      double d;
    };
    U u = {UINT64_C(0x3FF) << 52 | x >> 12};
    return u.d - 1.0;
  };

  uint64_t s[2];
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//! Class for maximization based on "Whale Optimization Algorithm"
//! doi:10.1016/j.advengsoft.2016.01.008
class EMMaximization
{
 private:
  EMRandom mRandom;
  size_t mPopSize;
  size_t mNIter;
  int mGenMethod;  // 0 = sky, 1 = cylinder, 2 = hspere
  double m_a;
  double m_a2;
  std::vector<std::vector<double> > mRanges;
  std::vector<std::vector<double> > mPopulation;
  double mBestCost;
  std::vector<double> mBestSolution;
  std::function<double(double, double)> mFunc;

 public:
  EMMaximization(const EMRandom& random, int genMethod)
    : mRandom(random)
    , mPopSize(200)
    , mNIter(500)
    , mGenMethod(genMethod)
    , m_a(0.)
    , m_a2(0.)
    , mBestCost(-1.)
  {
    mFunc = &DefaultJ;
  };
  ///////////////////////////////////////////////////////////////

  static double DefaultJ(double p, double theta)
  {
    double n = std::max(0.1, 2.856 - 0.655 * log(p));
    return 1600 * pow(p, 0.279) * pow(cos(theta), n);
  }
  ///////////////////////////////////////////////////////////////

  void SetParameters(double minP, double maxP, double minTheta, double maxTheta)
  {
    mRanges.push_back({minP, maxP});
    mRanges.push_back({minTheta, maxTheta});
  }
  ///////////////////////////////////////////////////////////////

  void SetParameters(double minP, double maxP, double minTheta, double maxTheta, double minPhi, double maxPhi)
  {
    mRanges.push_back({minP, maxP});
    mRanges.push_back({minTheta, maxTheta});
    mRanges.push_back({minPhi, maxPhi});
    mRanges.push_back({0, M_PI / 2.});
  }
  ///////////////////////////////////////////////////////////////

  void SetFunction(std::function<double(double, double)> func)
  {
    mFunc = func;
  }
  ///////////////////////////////////////////////////////////////

  double SkyFunc(double p, double theta)
  {
    return mFunc(p, theta) * cos(theta) * sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double CylFunc(double p, double theta)
  {
    return mFunc(p, theta) * pow(sin(theta), 2);
  }
  ///////////////////////////////////////////////////////////////

  double HSFunc(double p, double theta, double phi, double theta0)
  {
    return mFunc(p, theta) * (sin(theta0) * sin(theta) * cos(phi) + cos(theta0) * cos(theta)) * sin(theta);
  }
  ///////////////////////////////////////////////////////////////

  double Evaluate(std::vector<double>& v)
  {
    if (mGenMethod == 0)
    {
      return SkyFunc(v[0], v[1]);
    }
    else if (mGenMethod == 1)
    {
      return CylFunc(v[0], v[1]);
    }
    else
    {
      return HSFunc(v[0], v[1], v[2], v[3]);
    }
    return -1;
  }
  ///////////////////////////////////////////////////////////////

  void Evaluate()
  {
    double value;
    for (size_t i = 0; i < mPopSize; ++i)
    {
      value = Evaluate(mPopulation[i]);
      if (value > mBestCost)
      {
        mBestCost = value;
        mBestSolution = mPopulation[i];
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  void Init()
  {
    size_t dim = mRanges.size();
    mPopulation.resize(mPopSize);
    for (size_t i = 0; i < mPopSize; ++i)
    {
      mPopulation[i].resize(dim);
      for (size_t j = 0; j < dim; ++j)
      {
        mPopulation[i][j] = mRandom.GenerateRandomDouble(mRanges[j][0], mRanges[j][1]);
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  void UpdateParameters(size_t t)
  {
    m_a = 2. - t * (2. / mNIter);
    m_a2 = -1. + t * ((-1.) / mNIter);
  }
  ///////////////////////////////////////////////////////////////

  void Move()
  {
    double r1, r2, A, C, b, l, rw, p, D_tmp, D_best, distance;
    std::vector<double> tmp;
    for (size_t i = 0; i < mPopulation.size(); ++i)
    {
      r1 = mRandom.GenerateRandomDouble();
      r2 = mRandom.GenerateRandomDouble();
      A = 2 * m_a * r1 - m_a;
      C = 2 * r2;
      b = 1.;
      l = (m_a2 - 1) * mRandom.GenerateRandomDouble() + 1;
      p = mRandom.GenerateRandomDouble();

      for (size_t j = 0; j < mPopulation[0].size(); ++j)
      {
        if (p < 0.5)
        {
          if (fabs(A) >= 1)
          {
            rw = floor(mRandom.GenerateRandomDouble() * mPopulation.size());
            tmp = mPopulation[rw];
            D_tmp = fabs(C * tmp[j] - mPopulation[i][j]);
            mPopulation[i][j] = tmp[j] - A * D_tmp;
          }
          else
          {
            D_best = fabs(C * mBestSolution[j] - mPopulation[i][j]);
            mPopulation[i][j] = mBestSolution[j] - A * D_best;
          }
        }
        else
        {
          distance = fabs(mBestSolution[j] - mPopulation[i][j]);
          mPopulation[i][j] = distance * exp(b * l) * cos(l * 2 * M_PI) + mBestSolution[j];
        }
        if (mPopulation[i][j] < mRanges[j][0]) mPopulation[i][j] = mRanges[j][0];
        if (mPopulation[i][j] > mRanges[j][1]) mPopulation[i][j] = mRanges[j][1];
      }
    }
  }
  ///////////////////////////////////////////////////////////////

  double Maximize()
  {
    Init();
    Evaluate();
    for (size_t iter = 1; iter < mNIter; ++iter)
    {
      UpdateParameters(iter);
      Move();
      Evaluate();
    }
    return mBestCost;
  }
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//! Class for the generation of cosmic muons
class EcoMug
{
  friend class EMRandom;

 public:
  /// Possible generation methods
  enum EMGeometry
  {
    Sky,       ///< generation from a plane (flat sky)
    Cylinder,  ///< generation from a cylinder
    HSphere    ///< generation from a half-sphere
  };

 private:
  EMGeometry mGenMethod;
  std::array<double, 3> mGenerationPosition;
  double mGenerationTheta;
  double mGenerationPhi;
  double mGenerationMomentum;
  double mMinimumMomentum;
  double mMaximumMomentum;
  double mMinimumTheta;
  double mMaximumTheta;
  double mMinimumPhi;
  double mMaximumPhi;
  int mCharge;
  double mCylinderMinPositionPhi;
  double mCylinderMaxPositionPhi;
  double mHSphereMinPositionPhi;
  double mHSphereMaxPositionPhi;
  double mHSphereMinPositionTheta;
  double mHSphereMaxPositionTheta;
  double mHSphereCosMinPositionTheta;
  double mHSphereCosMaxPositionTheta;
  double mJPrime;
  double mN;
  double mRandAccRej;
  double mPhi0;
  double mTheta0;
  bool mAccepted;
  std::array<double, 2> mSkySize;
  std::array<double, 3> mSkyCenterPosition;
  double mCylinderHeight;
  double mCylinderRadius;
  std::array<double, 3> mCylinderCenterPosition;
  double mHSphereRadius;
  //  double mMaxFuncSkyCylinder;
  std::array<double, 3> mHSphereCenterPosition;
  EMRandom mRandom;
  std::array<double, 3> mMaxJ;
  std::array<double, 3> mMaxCustomJ;
  std::function<double(double, double)> mJ;

 public:
  // Default constructor
  EcoMug()
    : mGenMethod(Sky)
    , mGenerationPosition({{0., 0., 0.}})
    , mGenerationTheta(0.)
    , mGenerationPhi(0.)
    , mGenerationMomentum(0.)
    , mMinimumMomentum(0.01)
    , mMaximumMomentum(1000.)
    , mMinimumTheta(0.)
    , mMaximumTheta(M_PI / 2.)
    , mMinimumPhi(0.)
    , mMaximumPhi(2. * M_PI)
    , mCharge(1)
    , mCylinderMinPositionPhi(0.)
    , mCylinderMaxPositionPhi(2. * M_PI)
    , mHSphereMinPositionPhi(0.)
    , mHSphereMaxPositionPhi(2. * M_PI)
    , mHSphereMinPositionTheta(0.)
    , mHSphereMaxPositionTheta(M_PI / 2.)
    , mHSphereCosMinPositionTheta(1.)
    , mHSphereCosMaxPositionTheta(0.)
    , mJPrime(0.)
    , mN(0.)
    , mRandAccRej(0.)
    , mPhi0(0.)
    , mTheta0(0.)
    , mAccepted(false)
    , mSkySize({{0., 0.}})
    , mSkyCenterPosition({{0., 0., 0.}})
    , mCylinderHeight(0.)
    , mCylinderRadius(0.)
    , mCylinderCenterPosition({{0., 0., 0.}})
    , mHSphereRadius(0.)
    ,
    /* mMaxFuncSkyCylinder(5.3176), */ mHSphereCenterPosition({{0., 0., 0.}})
  {
    mMaxJ = {-1., -1., -1.};
    mMaxCustomJ = {-1., -1., -1.};
  };

  ///////////////////////////////////////////////////////////////
  // Methods to access the parameters of the generated muon
  ///////////////////////////////////////////////////////////////
  /// Get the generation position
  const std::array<double, 3>& GetGenerationPosition() const
  {
    return mGenerationPosition;
  };
  /// Get the generation momentum
  double GetGenerationMomentum() const
  {
    return mGenerationMomentum;
  };
  /// Get the generation momentum
  void GetGenerationMomentum(std::array<double, 3>& momentum) const
  {
    momentum = {
        mGenerationMomentum * sin(mGenerationTheta) * cos(mGenerationPhi),
        mGenerationMomentum * sin(mGenerationTheta) * sin(mGenerationPhi),
        mGenerationMomentum * cos(mGenerationTheta)};
  };
  /// Get the generation theta
  double GetGenerationTheta() const
  {
    return mGenerationTheta;
  };
  /// Get the generation phi
  double GetGenerationPhi() const
  {
    return mGenerationPhi;
  };
  /// Get charge
  int GetCharge() const
  {
    return mCharge;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Methods for the geometry of the generation
  ///////////////////////////////////////////////////////////////
  /// Set generation from sky
  void SetUseSky()
  {
    mGenMethod = Sky;
  };
  /// Set cylindrical generation
  void SetUseCylinder()
  {
    mGenMethod = Cylinder;
  };
  /// Set half-sphere generation
  void SetUseHSphere()
  {
    mGenMethod = HSphere;
  };
  /// Set the generation method (Sky, Cylinder or HSphere)
  void SetGenerationMethod(EMGeometry genM)
  {
    mGenMethod = genM;
  };

  /// Get the generation method (Sky, Cylinder or HSphere)
  EMGeometry GetGenerationMethod() const
  {
    return mGenMethod;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Common methods to all geometries
  ///////////////////////////////////////////////////////////////
  /// Set the differential flux J. Accepted functions are like
  /// double J(double momentum, double theta)
  /// momentum has to be in GeV/c and theta in radians
  void SetDifferentialFlux(std::function<double(double, double)> J)
  {
    mJ = J;
  };
  /// Set the seed for the internal PRNG (if 0 a random seed is used)
  void SetSeed(uint64_t seed)
  {
    if (seed > 0) mRandom.SetSeed(seed);
  };
  /// Set minimum generation Momentum
  void SetMinimumMomentum(double momentum)
  {
    mMinimumMomentum = momentum;
  };
  /// Set maximum generation Momentum
  void SetMaximumMomentum(double momentum)
  {
    mMaximumMomentum = momentum;
  };
  /// Set minimum generation Theta
  void SetMinimumTheta(double theta)
  {
    mMinimumTheta = theta;
  };
  /// Set maximum generation Theta
  void SetMaximumTheta(double theta)
  {
    mMaximumTheta = theta;
  };
  /// Set minimum generation Phi
  void SetMinimumPhi(double phi)
  {
    mMinimumPhi = phi;
  };
  /// Set maximum generation Phi
  void SetMaximumPhi(double phi)
  {
    mMaximumPhi = phi;
  };

  /// Get minimum generation Momentum
  double GetMinimumMomentum() const
  {
    return mMinimumMomentum;
  };
  /// Get maximum generation Momentum
  double GetMaximumMomentum() const
  {
    return mMaximumMomentum;
  };
  /// Get minimum generation Theta
  double GetMinimumTheta() const
  {
    return mMinimumTheta;
  };
  /// Get maximum generation Theta
  double GetMaximumTheta() const
  {
    return mMaximumTheta;
  };
  /// Get minimum generation Phi
  double GetMinimumPhi() const
  {
    return mMinimumPhi;
  };
  /// Get maximum generation Phi
  double GetMaximumPhi() const
  {
    return mMaximumPhi;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Methods for the plane-based generation
  ///////////////////////////////////////////////////////////////
  /// Set sky size
  void SetSkySize(const std::array<double, 2>& size)
  {
    mSkySize = size;
  };

  /// Set sky center position
  void SetSkyCenterPosition(const std::array<double, 3>& position)
  {
    mSkyCenterPosition = position;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Methods for the cylinder-based generation
  ///////////////////////////////////////////////////////////////
  /// Set cylinder radius
  void SetCylinderRadius(double radius)
  {
    mCylinderRadius = radius;
  };
  /// Set cylinder height
  void SetCylinderHeight(double height)
  {
    mCylinderHeight = height;
  };
  /// Set cylinder center position
  void SetCylinderCenterPosition(const std::array<double, 3>& position)
  {
    mCylinderCenterPosition = position;
  };
  void SetCylinderMinPositionPhi(double phi)
  {
    mCylinderMinPositionPhi = phi;
  };
  void SetCylinderMaxPositionPhi(double phi)
  {
    mCylinderMaxPositionPhi = phi;
  };
  /// Get cylinder radius
  double GetCylinderRadius() const
  {
    return mCylinderRadius;
  };
  /// Get cylinder height
  double GetCylinderHeight() const
  {
    return mCylinderHeight;
  };
  /// Get cylinder center position
  const std::array<double, 3>& GetCylinderCenterPosition() const
  {
    return mCylinderCenterPosition;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // Methods for the half sphere-based generation
  ///////////////////////////////////////////////////////////////
  /// Set half-sphere radius
  void SetHSphereRadius(double radius)
  {
    mHSphereRadius = radius;
  };
  /// Set half-sphere center position
  void SetHSphereCenterPosition(const std::array<double, 3>& position)
  {
    mHSphereCenterPosition = position;
  };
  void SetHSphereMinPositionPhi(double phi)
  {
    mHSphereMinPositionPhi = phi;
  };
  void SetHSphereMaxPositionPhi(double phi)
  {
    mHSphereMaxPositionPhi = phi;
  };
  void SetHSphereMinPositionTheta(double theta)
  {
    mHSphereMinPositionTheta = theta;
    mHSphereCosMinPositionTheta = cos(mHSphereMinPositionTheta);
  };
  void SetHSphereMaxPositionTheta(double theta)
  {
    mHSphereMaxPositionTheta = theta;
    mHSphereCosMaxPositionTheta = cos(mHSphereMaxPositionTheta);
  };
  /// Get half-sphere radius
  double GetHSphereRadius() const
  {
    return mHSphereRadius;
  };
  /// Get half-sphere center position
  const std::array<double, 3>& GetHSphereCenterPosition() const
  {
    return mHSphereCenterPosition;
  };
  ///////////////////////////////////////////////////////////////

 private:
  double F1Cumulative(double x)
  {
    return 1. - 8.534790171171021 / pow(x + 2.68, 87. / 40.);
  };

  double F1Inverse(double x)
  {
    return (2.68 - 2.68 * pow(1. - x, 40. / 87.)) / pow(1. - x, 40. / 87.);
  };

  double maxSkyJFunc()
  {
    return 1600 * pow(mMaximumMomentum, 0.279) * pow(cos(0.76158), 1.1) * sin(0.76158);
  };

  double maxCylJFunc()
  {
    return 1600 * pow(mMaximumMomentum, 0.279) * pow(cos(1.35081), 0.1) * pow(sin(1.35081), 2);
  };

  double maxHSJFunc()
  {
    return 1600 * pow(mMaximumMomentum, 0.279) * pow(cos(1.26452), 0.1) * (sin(1.26452) * sin(1.26452) + cos(1.26452) * cos(1.26452)) * sin(1.26452);
  };

  double GenerateMomentumF1()
  {
    double z = mRandom.GenerateRandomDouble(F1Cumulative(mMinimumMomentum), F1Cumulative(mMaximumMomentum));
    return F1Inverse(z);
  };

  void GeneratePositionSky()
  {
    mGenerationPosition[0] = mRandom.GenerateRandomDouble(mSkyCenterPosition[0] - mSkySize[0] / 2., mSkyCenterPosition[0] + mSkySize[0] / 2.);
    mGenerationPosition[1] = mRandom.GenerateRandomDouble(mSkyCenterPosition[1] - mSkySize[1] / 2., mSkyCenterPosition[1] + mSkySize[1] / 2.);
    mGenerationPosition[2] = mSkyCenterPosition[2];
  };

  void GeneratePositionCylinder()
  {
    mPhi0 = mRandom.GenerateRandomDouble(mCylinderMinPositionPhi, mCylinderMaxPositionPhi);
    mGenerationPosition[0] = mCylinderCenterPosition[0] + mCylinderRadius * cos(mPhi0);
    mGenerationPosition[1] = mCylinderCenterPosition[1] + mCylinderRadius * sin(mPhi0);
    mGenerationPosition[2] = mRandom.GenerateRandomDouble(mCylinderCenterPosition[2] - mCylinderHeight / 2., mCylinderCenterPosition[2] + mCylinderHeight / 2.);
  };

  void ComputeMaximumCustomJ()
  {
    EMMaximization maximizer(mRandom, mGenMethod);
    maximizer.SetFunction(mJ);
    if (mGenMethod == 0 || mGenMethod == 1)
    {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta);
    }
    else
    {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta, mMinimumPhi, mMaximumPhi);
    }
    mMaxCustomJ[mGenMethod] = maximizer.Maximize();
  };

  void ComputeMaximum()
  {
    EMMaximization maximizer(mRandom, mGenMethod);
    if (mGenMethod == 0 || mGenMethod == 1)
    {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta);
    }
    else
    {
      maximizer.SetParameters(mMinimumMomentum, mMaximumMomentum, mMinimumTheta, mMaximumTheta, mMinimumPhi, mMaximumPhi);
    }
    mMaxJ[mGenMethod] = maximizer.Maximize();
  };

 public:
  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon from the pre-defined J
  ///////////////////////////////////////////////////////////////
  void Generate()
  {
    mAccepted = false;

    if (mMaxJ[mGenMethod] < 0) ComputeMaximum();

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder)
    {
      // Generation of the momentum and theta angle
      while (!mAccepted)
      {
        mRandAccRej = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = GenerateMomentumF1();
        mN = 2.856 - 0.655 * log(mGenerationMomentum);
        if (mN < 0.1) mN = 0.1;

        if (mGenMethod == Sky)
        {
          mJPrime = 1600 * pow(mGenerationMomentum, 0.279) * pow(cos(mGenerationTheta), mN + 1) * sin(mGenerationTheta);
          if (mMaxJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
        }

        if (mGenMethod == Cylinder)
        {
          mJPrime = 1600 * pow(mGenerationMomentum, 0.279) * pow(cos(mGenerationTheta), mN) * pow(sin(mGenerationTheta), 2);
          if (mMaxJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky)
      {
        GeneratePositionSky();
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      }
      if (mGenMethod == Cylinder)
      {
        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted)
        {
          mRandAccRej = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2. * M_PI) mGenerationPhi -= 2. * M_PI;

        // Check if the muon is inward
        if (sin(mGenerationTheta) * cos(mGenerationPhi) * mGenerationPosition[0] + sin(mGenerationTheta) * sin(mGenerationPhi) * mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere)
    {
      // Generation point on the half-sphere
      mPhi0 = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted)
      {
        mRandAccRej = mRandom.GenerateRandomDouble();
        mTheta0 = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = GenerateMomentumF1();
        mN = 2.856 - 0.655 * log(mGenerationMomentum);
        if (mN < 0.1) mN = 0.1;

        mJPrime = 1600 * pow(mGenerationMomentum, 0.279) * pow(cos(mGenerationTheta), mN) * (sin(mGenerationTheta) * sin(mTheta0) * cos(mGenerationPhi) + cos(mGenerationTheta) * cos(mTheta0)) * sin(mGenerationTheta);
        if (mMaxJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
      }

      mGenerationPosition[0] = mHSphereRadius * sin(mTheta0) * cos(mPhi0) + mHSphereCenterPosition[0];
      mGenerationPosition[1] = mHSphereRadius * sin(mTheta0) * sin(mPhi0) + mHSphereCenterPosition[1];
      mGenerationPosition[2] = mHSphereRadius * cos(mTheta0) + mHSphereCenterPosition[2];

      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + mPhi0;
      if (mGenerationPhi >= 2 * M_PI) mGenerationPhi -= 2 * M_PI;

      mGenerationPhi += M_PI;
      if (mGenerationPhi >= 2 * M_PI) mGenerationPhi -= 2 * M_PI;
    }

    // Generate the charge
    if (mRandom.GenerateRandomDouble(-100,128)>=0) mCharge = 1;
    else mCharge = -1;
  };
  ///////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  /// Generate a cosmic muon for the user-defined J
  ///////////////////////////////////////////////////////////////
  void GenerateFromCustomJ()
  {
    mAccepted = false;

    if (mMaxCustomJ[mGenMethod] < 0) ComputeMaximumCustomJ();

    // Sky or cylinder generation
    if (mGenMethod == Sky || mGenMethod == Cylinder)
    {
      // Generation of the momentum and theta angle
      while (!mAccepted)
      {
        mRandAccRej = mRandom.GenerateRandomDouble();
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);

        if (mGenMethod == Sky)
        {
          mJPrime = mJ(mGenerationMomentum, mGenerationTheta) * cos(mGenerationTheta) * sin(mGenerationTheta);
          if (mMaxCustomJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
        }

        if (mGenMethod == Cylinder)
        {
          mJPrime = mJ(mGenerationMomentum, mGenerationTheta) * pow(sin(mGenerationTheta), 2) * cos(mGenerationPhi);
          if (mMaxCustomJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
        }
      }
      mGenerationTheta = M_PI - mGenerationTheta;

      // Generation of the position and phi angle
      if (mGenMethod == Sky)
      {
        GeneratePositionSky();
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
      }
      if (mGenMethod == Cylinder)
      {
        mAccepted = false;
        GeneratePositionCylinder();
        while (!mAccepted)
        {
          mRandAccRej = mRandom.GenerateRandomDouble();
          mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
          if (mRandAccRej < fabs(cos(mGenerationPhi))) mAccepted = true;
        }
        mGenerationPhi = mGenerationPhi + mPhi0;
        if (mGenerationPhi >= 2. * M_PI) mGenerationPhi -= 2. * M_PI;

        // Check if the muon is inward
        if (sin(mGenerationTheta) * cos(mGenerationPhi) * mGenerationPosition[0] + sin(mGenerationTheta) * sin(mGenerationPhi) * mGenerationPosition[1] > 0) Generate();
      }
    }

    // Half-sphere generation
    if (mGenMethod == HSphere)
    {
      // Generation point on the half-sphere
      mPhi0 = mRandom.GenerateRandomDouble(mHSphereMinPositionPhi, mHSphereMaxPositionPhi);
      while (!mAccepted)
      {
        mRandAccRej = mRandom.GenerateRandomDouble();
        mTheta0 = acos(mRandom.GenerateRandomDouble(mHSphereCosMaxPositionTheta, mHSphereCosMinPositionTheta));
        mGenerationTheta = mRandom.GenerateRandomDouble(mMinimumTheta, mMaximumTheta);
        mGenerationPhi = mRandom.GenerateRandomDouble(mMinimumPhi, mMaximumPhi);
        mGenerationMomentum = mRandom.GenerateRandomDouble(mMinimumMomentum, mMaximumMomentum);

        mJPrime = mJ(mGenerationMomentum, mGenerationTheta) * (sin(mTheta0) * sin(mGenerationTheta) * cos(mGenerationPhi) + cos(mTheta0) * cos(mGenerationTheta)) * sin(mGenerationTheta);
        if (mMaxCustomJ[mGenMethod] * mRandAccRej < mJPrime) mAccepted = true;
      }

      mGenerationPosition[0] = mHSphereRadius * sin(mTheta0) * cos(mPhi0) + mHSphereCenterPosition[0];
      mGenerationPosition[1] = mHSphereRadius * sin(mTheta0) * sin(mPhi0) + mHSphereCenterPosition[1];
      mGenerationPosition[2] = mHSphereRadius * cos(mTheta0) + mHSphereCenterPosition[2];

      mGenerationTheta = M_PI - mGenerationTheta;
      mGenerationPhi = mGenerationPhi + mPhi0;
      if (mGenerationPhi >= 2 * M_PI) mGenerationPhi -= 2 * M_PI;

      mGenerationPhi += M_PI;
      if (mGenerationPhi >= 2 * M_PI) mGenerationPhi -= 2 * M_PI;
    }

    // Generate the charge
    if (mRandom.GenerateRandomDouble(-100,128)>=0) mCharge = 1;
    else mCharge = -1;
  };
  ///////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#endif
