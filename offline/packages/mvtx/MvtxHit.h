/**
 * @file mvtx/MvtxHit.h
 * @author D. McGlinchey
 * @date June 2018
 * @brief Mvtx hit object
 */
#ifndef MVTX_MVTXHIT_H
#define MVTX_MVTXHIT_H

#include <trackbase/TrkrHit.h>

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif

/**
 * @brief Mvtx hit object
 *
 * Container for Mvtx hit object representing
 * a single hit pixel within a chip
 */
class MvtxHit : public TrkrHit
{
 public:
  //! ctor
  MvtxHit();
  //! dtor
  virtual ~MvtxHit(){};

  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual int isValid() const;

 private:
  ClassDef(MvtxHit, 1);
};

#endif  //MVTX_MVTXHIT_H
