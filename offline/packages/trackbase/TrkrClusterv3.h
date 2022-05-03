/**
 * @file trackbase/TrkrClusterv3.h
 * @author J. Osborn
 * @date October 2021
 * @brief Version 3 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV3_H
#define TRACKBASE_TRKRCLUSTERV3_H

#include "TrkrCluster.h"
#include "TrkrDefs.h"
#include <iostream>

class PHObject;

/**
 * @brief Version 3 of TrkrCluster
 *
 * This version of TrkrCluster contains Acts source link objects
 * as member variables, to join the Svtx and Acts worlds
 */
class TrkrClusterv3 : public TrkrCluster
{
 public:
 
  //! ctor
  TrkrClusterv3();

  //!dtor
  ~TrkrClusterv3() override = default;
  
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new TrkrClusterv3(*this); }
  
  //! copy content from base class
  void CopyFrom( const TrkrCluster& ) override;

  //! copy content from base class
  void CopyFrom( TrkrCluster* source ) override
  { CopyFrom( *source ); } 
  
  //
  // cluster position
  //
  float getPosition(int coor) const override { return m_local[coor]; }
  void setPosition(int coor, float xi) override { m_local[coor] = xi; }
  float getLocalX() const override { return m_local[0]; }
  void setLocalX(float loc0) override { m_local[0] = loc0; }
  float getLocalY() const override { return m_local[1]; }
  void setLocalY(float loc1) override { m_local[1] = loc1; }

  /// Acts functions, for Acts module use only
  void setActsLocalError(unsigned int i, unsigned int j, float value) override;
  float getActsLocalError(unsigned int i, unsigned int j) const override { return m_actsLocalErr[i][j]; }
  TrkrDefs::subsurfkey getSubSurfKey() const override { return m_subsurfkey; }
  void setSubSurfKey(TrkrDefs::subsurfkey id) override { m_subsurfkey = id; }

  //
  // cluster info
  //
  unsigned int getAdc() const override { return m_adc; }
  void setAdc(unsigned int adc) override { m_adc = adc; }

  //
  // convenience interface
  //
  float getRPhiError() const override;
  float getZError() const override;

  /// deprecated global funtions with a warning
  float getX() const override 
  { std::cout << "Deprecated getx trkrcluster function!"<<std::endl; return NAN;}
  float getY() const override 
  { std::cout << "Deprecated gety trkrcluster function!"<<std::endl; return NAN;}
  float getZ() const override 
  { std::cout << "Deprecated getz trkrcluster function!"<<std::endl; return NAN;}
   void setX(float) override 
   { std::cout << "Deprecated setx trkrcluster function!"<<std::endl;} 
   void setY(float) override 
   { std::cout << "Deprecated sety trkrcluster function!"<<std::endl;} 
   void setZ(float) override 
   { std::cout << "Deprecated setz trkrcluster function!"<<std::endl;}
   float getSize(unsigned int, unsigned int) const override 
   {std::cout << "Deprecated getsize trkrcluster function!" << std::endl; return NAN;}       
   void setSize(unsigned int, unsigned int, float) override 
   {std::cout << "Deprecated setsize trkrcluster function!" << std::endl;}
   float getError(unsigned int, unsigned int) const override 
   {std::cout << "Deprecated geterr trkrcluster function!" << std::endl; return NAN;}
   void setError(unsigned int, unsigned int, float) override 
   { std::cout << "Deprecated seterr trkrcluster function!" << std::endl; }
   float getPhiSize() const override 
   { std::cout << "Deprecated size function"<< std::endl; return NAN;}
   float getZSize() const override 
   {std::cout << "Deprecated size function" << std::endl; return NAN;}
   float getPhiError() const override 
   { std::cout << "Deprecated getPhiError function"<< std::endl; return NAN;}

 protected:

  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  TrkrDefs::subsurfkey m_subsurfkey; //< unique identifier for hitsetkey-surface maps

  unsigned int m_adc;           //< cluster sum adc (D. McGlinchey - Do we need this?)
  
  float m_local[2];             //< 2D local position [cm]
  float m_actsLocalErr[2][2];   //< 2D local error for Acts [cm]

  ClassDefOverride(TrkrClusterv3, 2)
};

#endif //TRACKBASE_TRKRCLUSTERV3_H
