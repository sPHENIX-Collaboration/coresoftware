#ifndef __BBCOUT_H__
#define __BBCOUT_H__

#include "phool/PHObject.h"

///
class BbcOut: public PHObject
{
public:
  ///
  virtual ~BbcOut() {}
  
  /** identify Function from PHObject
      @param os Output Stream 
   */
  virtual void identify(std::ostream& os = std::cout) const override; 

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const override;

  /// get ZVertex determined by Bbc
  virtual Float_t get_VertexPoint() const;

  /// get Error on ZVertex determined by Bbc
  virtual Float_t get_dVertexPoint() const;

  /// get T0 determined by Bbc
  virtual Float_t get_TimeZero() const;

  /// get Error on T0 determined by Bbc
  virtual Float_t get_dTimeZero() const;

  /** set T0, Error on T0, ZVertex and Error on ZVertex
      @param t0 Bbc T0
      @param t0err Bbc Error on T0
      @param vtx Bbc ZVertex 
      @param vtxerr Bbc Error on ZVertex
   */
  virtual void set_TimeVertex(const Float_t t0, const Float_t t0err, const Float_t vtx, const Float_t vtxerr)
  { 
    set_TimeZero( t0, t0err );
    set_Vertex( vtx, vtxerr );
  }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  virtual void set_TimeZero(const Float_t t0, const Float_t t0err = 0);

  //! set vertex
  virtual void set_Vertex( const Float_t vtx, const Float_t vtxerr);
  
  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
   */
  virtual void set_dZVertex(const Float_t vtxerr);
  
  /** Add Bbc North/South object containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param ncharge Number of Charged Particles into North/South
      @param timing Timing of North/South
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual void AddBbcNS(const int iBBC, const Short_t npmt, const Float_t ncharge, const Float_t timing);

  /** get Number of PMT's fired in North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual Short_t get_nPMT(const int iBBC) const;

  /** get Number of Charged Particles into North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual Float_t get_nCharge(const int iBBC) const;

  /** get Timing of North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual Float_t get_Timing(const int iBBC) const;

  virtual void FillFromClass(const BbcOut& old);
    
private:
  void virtual_warning(const char *funcname) const;

  /// Root Internal Version
  ClassDefOverride(BbcOut,1)


};

#endif

