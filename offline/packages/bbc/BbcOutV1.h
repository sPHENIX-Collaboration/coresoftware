// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUTV1_H
#define BBC_BBCOUTV1_H

#include "BbcOut.h"

class TClonesArray;

///
class BbcOutV1: public BbcOut
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
  void identify(std::ostream& os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  Float_t get_VertexPoint() const override {return Bbc_ZVertex;}

  /// get Error on ZVertex determined by Bbc
  Float_t get_dVertexPoint() const override {return Bbc_dZVertex;}

  /// get T0 determined by Bbc
  Float_t get_TimeZero() const override {return Bbc_TimeZero;}

  /// get Error on T0 determined by Bbc
  Float_t get_dTimeZero() const override {return Bbc_dTimeZero;}

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  void set_TimeZero(const Float_t t0, const Float_t t0err = 0) override;

  //! set vertex
  void set_Vertex( const Float_t vtx, const Float_t vtxerr = 0) override;
  
  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
  */
  void set_dZVertex(const Float_t vtxerr) override;

  /** Add Bbc North/South object containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param energy Energy in North/South
      @param timing Timing of North/South
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  void AddBbcNS(const int nBbc, const Short_t npmt, const Float_t chargesum, const Float_t timing) override;

  /** get Number of PMT's fired in North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  Short_t get_nPMT(const int nBbc) const override;


  /** get Number of Charged Particles into North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  Float_t get_nCharge(const int nBbc) const override;

  /** get Timing of North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  Float_t get_Timing(const int nBbc) const override;

  virtual void Init();

private:
  
  
  TClonesArray *GetBbcNS() const { return BbcNS; }

  Float_t Bbc_ZVertex;
  Float_t Bbc_dZVertex;
  Float_t Bbc_TimeZero;
  Float_t Bbc_dTimeZero;
  TClonesArray *BbcNS;

  ClassDefOverride(BbcOutV1,1)

};

#endif
