<<<<<<< HEAD
#ifndef MINIMUMBIASINFO_H
#define MINIMUMBIASINFO_H
=======
#ifndef TRIGGER_MINIMUMBIASINFO_H
#define TRIGGER_MINIMUMBIASINFO_H
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

#include <phool/PHObject.h>

class MinimumBiasInfo : public PHObject
{
 public:
  ~MinimumBiasInfo() override{};

  void identify(std::ostream &os = std::cout) const override { os << "MinimumBiasInfo base class" << std::endl; };
<<<<<<< HEAD
  virtual void Reset() override {}
  int isValid() const override { return 0; }

 protected:

=======
  void Reset() override {}
  int isValid() const override { return 0; }
  virtual void setIsAuAuMinimumBias(bool) { return; }
  virtual bool isAuAuMinimumBias() const { return false; }

 protected:
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  MinimumBiasInfo() {}

 private:
  ClassDefOverride(MinimumBiasInfo, 1);
};

<<<<<<< HEAD
#endif  // TRIGGER_MINBIASTRIGGERINFO_H
=======
#endif  // TRIGGER_MINIMUMBIASINFO_H
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
