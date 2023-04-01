// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUTV1_H
#define BBC_BBCOUTV1_H

#include "BbcOut.h"

#include <iostream>

class TClonesArray;

///
class BbcOutV1 : public BbcOut
{
 public:
  ///
  BbcOutV1();
  ///
  ~BbcOutV1() override;

  /// Clear Event from memory
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  float get_VertexPoint() const override { return Bbc_ZVertex; }

  /// get Error on ZVertex determined by Bbc
  float get_dVertexPoint() const override { return Bbc_dZVertex; }

  /// get T0 determined by Bbc
  float get_TimeZero() const override { return Bbc_TimeZero; }

  /// get Error on T0 determined by Bbc
  float get_dTimeZero() const override { return Bbc_dTimeZero; }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  void set_TimeZero(const float t0, const float t0err = 0) override;

  //! set vertex
  void set_Vertex(const float vtx, const float vtxerr = 0) override;

  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
  */
  void set_dZVertex(const float vtxerr) override;

  /** Add Bbc North/South object containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param energy Energy in North/South
      @param timing Timing of North/South
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  void AddBbcNS(const int nBbc, const short npmt, const float chargesum, const float timing) override;

  /** get Number of PMT's fired in North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  short get_nPMT(const int nBbc) const override;

  /** get Number of Charged Particles into North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  float get_nCharge(const int nBbc) const override;

  /** get Timing of North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  float get_Timing(const int nBbc) const override;


 private:

  TClonesArray *GetBbcNS() const { return BbcNS; }
  void Init();

  float Bbc_ZVertex{};
  float Bbc_dZVertex{};
  float Bbc_TimeZero{};
  float Bbc_dTimeZero{};
  TClonesArray *BbcNS;

  ClassDefOverride(BbcOutV1, 1)
};

#endif
