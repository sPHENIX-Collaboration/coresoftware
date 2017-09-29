#ifndef __TPCELECTRONPDF_H__
#define __TPCELECTRONPDF_H__

#include <vector>
#include <cmath>
#include <TPCHit.h>

class TPCElectronPDF : public TPCHit {
 public:
  TPCElectronPDF();
  virtual ~TPCElectronPDF();

  void Reset();
  unsigned int GetN() {return fW.size();}
  void AddElectron(float w);
  float GetElectrons(unsigned int i) {return i<GetN()?fW[i]:0;}
  void Amplify(unsigned int i, float v) {if(i<GetN()) fW[i]*=v;}
  float GetRMST(unsigned int i) {return i<GetN()?std::sqrt(fMST[i]):0;}
  float GetRMSL0(unsigned int i) {return i<GetN()?std::sqrt(fMSL0[i]):0;}
  float GetRMSL1(unsigned int i) {return i<GetN()?std::sqrt(fMSL1[i]):0;}
  void AddMST(unsigned int i, float v) {if(i<GetN()) fMST[i]+=v;}
  void AddMSL0(unsigned int i, float v) {if(i<GetN()) fMSL0[i]+=v;}
  void AddMSL1(unsigned int i, float v) {if(i<GetN()) fMSL1[i]+=v;}
  void SetMSL0(unsigned int i, float v) {if(i<GetN()) fMSL0[i]=v;}
  void SetMSL1(unsigned int i, float v) {if(i<GetN()) fMSL1[i]=v;}
  virtual void CopyFrom(TPCHit *hit);
  float GetR(unsigned int i);
  float GetPhi(unsigned int i);
  float GetZ() {return TPCHit::GetZ();}
  float GetZ(unsigned int i) {return i<GetN()?TPCHit::GetZ()+fDZ[i]:TPCHit::GetZ();}
  float GetT(unsigned int i) {return i<GetN()?fT[i]:0;}
  float GetDx(unsigned int i) {return i<GetN()?fDX[i]:0;}
  float GetDy(unsigned int i) {return i<GetN()?fDY[i]:0;}
  float GetDz(unsigned int i) {return i<GetN()?fDZ[i]:0;}
  void SetT(unsigned int i, float v) {if(i<GetN()) fT[i] = v;}
  void AddDx(unsigned int i, float v) {if(i<GetN()) fDX[i] += v;}
  void AddDy(unsigned int i, float v) {if(i<GetN()) fDY[i] += v;}
  void AddDz(unsigned int i, float v) {if(i<GetN()) fDZ[i] += v;}

 protected:
  std::vector<float> fW; // weight
  std::vector<float> fDX; // distortions in Rad
  std::vector<float> fDY; // distortions in Phi
  std::vector<float> fDZ; // distortions in Z
  std::vector<float> fT;  // stores times
  std::vector<float> fMST; // mean square in tranverse direction
  std::vector<float> fMSL0; // mean square in longudinal direction
  std::vector<float> fMSL1; // mean square in longudinal direction
};

#endif
