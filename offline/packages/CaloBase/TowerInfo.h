#ifndef TowerInfo_H
#define TowerInfo_H

#include <cstddef>  // for size_t
#include <iostream>
#include <set>
#include <TObject.h>

#include <phool/PHObject.h>

class TowerInfo : public TObject
{
 public:
  TowerInfo();
  TowerInfo(const TowerInfo &ti);
  ~TowerInfo() override;

  void setTime(short t) { _time = t; }
  float getTime() { return _time; }
  void setAmplitude(float amp) { _amplitude = amp; }
  float getAmplitude() { return _amplitude; }


protected:
  short _time;
  float _amplitude;


  ClassDefOverride(TowerInfo, 1);
};

#endif
