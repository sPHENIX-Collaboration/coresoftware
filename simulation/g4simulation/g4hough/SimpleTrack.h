#ifndef __SIMPLETRACK__
#define __SIMPLETRACK__

#include "TObject.h"
#include <vector>


class SimpleMCHit : public TObject
{
public:
  SimpleMCHit(float xx=0., float yy=0., float zz=0., unsigned int ll=0, unsigned int ii=0) : x(xx), y(yy), z(zz), layer(ll), index(ii) {}
  ~SimpleMCHit(){}
  
  float x;
  float y;
  float z;
  unsigned int layer;
  unsigned int index;
  
  ClassDef(SimpleMCHit, 1);
};


class SimpleMCTrack : public TObject
{
  public:
    SimpleMCTrack() : px(0.), py(0.), pz(0.), d(0.), z0(0.)  {}
    ~SimpleMCTrack(){}
    
    std::vector<SimpleMCHit> hits;
    float px;
    float py;
    float pz;
    float d;
    float z0;
    
    ClassDef(SimpleMCTrack, 1);
};


class SimpleMCEvent : public TObject
{
  public:
    SimpleMCEvent(){}
    ~SimpleMCEvent(){}
    
    std::vector<SimpleMCTrack> tracks;
    
    ClassDef(SimpleMCEvent, 1);
};


class SimpleRecoTrack : public TObject
{
  public:
    SimpleRecoTrack() : px(0.), py(0.), pz(0.), d(0.), z0(0.), quality(0.), charge(0.) {}
    ~SimpleRecoTrack(){}
    
    std::vector<unsigned int> indexes;
    float px;
    float py;
    float pz;
    float d;
    float z0;
    float quality;
    short int charge;
    
    ClassDef(SimpleRecoTrack, 1);
};


class SimpleRecoEvent : public TObject
{
  public:
    SimpleRecoEvent(){}
    ~SimpleRecoEvent(){}
    
    std::vector<SimpleRecoTrack> tracks;
    
    ClassDef(SimpleRecoEvent, 1);
};

#endif

