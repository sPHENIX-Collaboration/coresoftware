#ifndef TRAININGHITSCONTAINER_H
#define TRAININGHITSCONTAINER_H

#include "TrainingHits.h"
#include <vector>

class TrainingHitsContainer: public PHObject
{
  public:
    TrainingHitsContainer();
    ~TrainingHitsContainer() override {}
    void Reset() override;

    std::vector<TrainingHits> v_hits;

    ClassDefOverride(TrainingHitsContainer, 1)
};

#endif //TRAININGHITSCONTAINER_H
