/**
 * @file trackbase/TrkrClusterv4.h
 * @author C. Roland
 * @date April 2022
 * @brief Version 4 of TrkrCluster
 */
#ifndef TRACKBASE_TRKRCLUSTERV4_H
#define TRACKBASE_TRKRCLUSTERV4_H

#include "TrkrCluster.h"
#include "TrkrDefs.h"
#include <iostream>

class PHObject;

/**
 * @brief Version 4 of TrkrCluster
 *
 * This version of TrkrCluster is reduced to a minimal number of data members
 */
class TrkrClusterv4 : public TrkrCluster
{
 public:
 
  //! ctor
  TrkrClusterv4();

  //!dtor
  ~TrkrClusterv4() override = default;
  
  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override {}
  int isValid() const override;
  PHObject* CloneMe() const override { return new TrkrClusterv4(*this); }
  
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
  float getRPhiError() const override
  { std::cout << "Deprecated getRPhiError trkrcluster function!"<<std::endl; return NAN;}
  float getZError() const override
  { std::cout << "Deprecated getZError trkrcluster function!"<<std::endl; return NAN;}

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

   char getSize() { return m_size; }
   void setSize(char size) { m_size = size; }
 
   float getPhiSize() const override { return (float) m_phisize; }
   void setPhiSize(char phisize) { m_phisize = phisize; }
 
   float getZSize() const override { return (float) m_zsize; }
   void setZSize(char zsize) { m_zsize = zsize; }
 
   char getOverlap() { return m_overlap; }
   void setOverlap(char overlap) { m_overlap = overlap; }
 
   char getEdge() { return m_edge; }
   void setEdge(char edge) { m_edge = edge; }

   //float getPhiSize() const override 
   //{ std::cout << "Deprecated size function"<< std::endl; return NAN;}
   //float getZSize() const override 
   //{std::cout << "Deprecated size function" << std::endl; return NAN;}
   //float getPhiError() const override 
   //{ std::cout << "Deprecated getPhiError function"<< std::endl; return NAN;}

 protected:

  TrkrDefs::cluskey m_cluskey;  //< unique identifier within container
  TrkrDefs::subsurfkey m_subsurfkey; //< unique identifier for hitsetkey-surface maps

  unsigned short int m_adc;           //< cluster sum adc
  
  float m_local[2];          //< 2D local position [cm]
  bool m_valid;
  char m_size;
  char m_phisize;
  char m_zsize;
  char m_overlap;
  char m_edge;

  ClassDefOverride(TrkrClusterv4, 2)
};

#endif //TRACKBASE_TRKRCLUSTERV4_H
