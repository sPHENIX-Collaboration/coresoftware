// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUT_H
#define BBC_BBCOUT_H

#include <phool/PHObject.h>

#include <string>

///
class BbcOut : public PHObject
{
 public:
  ///
  ~BbcOut() override {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  void Reset() override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  virtual float get_VertexPoint() const;

  /// get Error on ZVertex determined by Bbc
  virtual float get_dVertexPoint() const;

  /// get T0 determined by Bbc
  virtual float get_TimeZero() const;

  /// get Error on T0 determined by Bbc
  virtual float get_dTimeZero() const;

  /** set T0, Error on T0, ZVertex and Error on ZVertex
      @param t0 Bbc T0
      @param t0err Bbc Error on T0
      @param vtx Bbc ZVertex
      @param vtxerr Bbc Error on ZVertex
   */
  virtual void set_TimeVertex(const float t0, const float t0err, const float vtx, const float vtxerr)
  {
    set_TimeZero(t0, t0err);
    set_Vertex(vtx, vtxerr);
  }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  virtual void set_TimeZero(const float t0, const float t0err = 0);

  //! set vertex
  virtual void set_Vertex(const float vtx, const float vtxerr);

  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
   */
  virtual void set_dZVertex(const float vtxerr);

  /** Add Bbc North/South object containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param ncharge Number of Charged Particles into North/South
      @param timing Timing of North/South
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual void AddBbcNS(const int iBBC, const short npmt, const float ncharge, const float timing);

  /** get Number of PMT's fired in North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual short get_nPMT(const int iBBC) const;

  /** get Number of Charged Particles into North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual float get_nCharge(const int iBBC) const;

  /** get Timing of North/South Bbc
      @param nBbc  Arm, use Bbc::North and Bbc::South
   */
  virtual float get_Timing(const int iBBC) const;

  virtual void FillFromClass(const BbcOut& old);

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(BbcOut, 1)
};

#endif
