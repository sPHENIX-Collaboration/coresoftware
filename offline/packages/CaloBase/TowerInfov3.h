#ifndef TOWERINFOV3_H
#define TOWERINFOV3_H

#include "TowerInfo.h"
#include "TowerInfov2.h"

#include <cstdint> // For int16_t

class TowerInfov3: public TowerInfov2
{
public:
    TowerInfov3() {}
    ~TowerInfov3() override {}
    
    void Reset() override;
    void Clear(Option_t* = "") override;
    
    // Getter and setter for waveform
    int get_nsample() const override {return nsample;}
    int16_t get_waveform_value(int index) const override;
    void set_waveform_value(int index, int16_t value) override;

    void copy_tower(TowerInfo* tower) override;
    
private:
    static const int nsample = 31;
    int16_t _waveform[nsample] = {0}; // Initializes the entire array to zero

    ClassDefOverride(TowerInfov3, 1);
    // Inherit other methods and properties from TowerInfov2
};

#endif
