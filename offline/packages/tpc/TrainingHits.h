#ifndef TRAININGHITS_H
#define TRAININGHITS_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <array>

class TrainingHits : public PHObject
{
 public:
  TrainingHits();
  ~TrainingHits() override {}
  void Reset() override;

  static const int nd = 5;
  std::array<Short_t, (2 * nd + 1) * (2 * nd + 1)> v_adc{};
  Float_t radius;
  Float_t phi;
  Float_t z;
  Float_t phistep;
  Float_t zstep;
  Short_t layer;
  Short_t ntouch;
  Short_t nedge;
  TrkrDefs::cluskey cluskey;

  ClassDefOverride(TrainingHits, 1)
};

#endif  // TRAININGHITS_H
