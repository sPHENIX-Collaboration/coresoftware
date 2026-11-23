#ifndef KFPARTICLESPHENIX_KFPARTICLEUNSORTED_H
#define KFPARTICLESPHENIX_KFPARTICLEUNSORTED_H

#include <string>

// TEMPORARY. Holds the code taken out of KFParticle_Tools that KFParticle_eventReconstruction does not use,
// while it is redistributed during the Tools rework. Nothing new should be added here, and this file
// should be gone by the end of the rework.
class KFParticle_unsorted
{
 public:
  KFParticle_unsorted() = delete;

  // Only used by KFParticle_sPHENIX::parseDecayDescriptor
  static bool findParticle(const std::string &particle);
};

#endif  // KFPARTICLESPHENIX_KFPARTICLEUNSORTED_H
