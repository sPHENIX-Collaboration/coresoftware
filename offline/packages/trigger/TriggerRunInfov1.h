#ifndef TRIGGER_TRIGGERRUNINFOV1_H_
#define TRIGGER_TRIGGERRUNINFOV1_H_

#include "TriggerRunInfo.h"

#include <array>
#include <cstdint>
#include <iostream>
#include <string>

class TriggerRunInfov1 : public TriggerRunInfo
{
 public:
  TriggerRunInfov1() = default;

  ~TriggerRunInfov1() override = default;

  void setTrigger(int index, const std::string& name, int bit, int prescale) override;

  void setTriggerScalers(int index, int scalertype, uint64_t scalers) override;

  void setTriggerPrescale(int index, double prescale) override;

  double getPrescaleByName(const std::string& name) const override;

  double getPrescaleByBit(int triggerbit) const override;

  int getInitialPrescaleByName(const std::string& name) const override;

  int getInitialPrescaleByBit(int triggerbit) const override;

  uint64_t getScalersByName(const std::string& name) const override;

  uint64_t getScalersByBit(int triggerbit) const override;

  uint64_t getLiveScalersByName(const std::string& name) const override;

  uint64_t getLiveScalersByBit(int triggerbit) const override;

  uint64_t getRawScalersByName(const std::string& name) const override;

  uint64_t getRawScalersByBit(int triggerbit) const override;

  uint32_t getTriggerBitByName(const std::string& name) const override;

  void identify(std::ostream& os = std::cout) const override;

  std::string getTriggerName(int triggerbit) const override;

 private:
  std::array<std::string, 64> trigger_names{};
  std::array<int, 64> trigger_bits{};
  std::array<int, 64> trigger_initial_prescales{};
  std::array<double, 64> trigger_prescales{};
  std::array<std::array<uint64_t, 3>, 64> trigger_scalers{};

 private:  // so the ClassDef does not show up with doc++
  ClassDefOverride(TriggerRunInfov1, 1);
};

#endif
