// @file mvtx_utils.h
// @brief Declarations of helper classes for the MVTX

#ifndef MVTXDECODER_UTILS_H
#define MVTXDECODER_UTILS_H

#include <cstdint>
#include <cassert>
#include <iostream>
#include <limits>
#include <array>
#include <memory>
#include <phool/recoConsts.h>

namespace mvtx_utils {

#define clean_errno() (errno == 0 ? "None" : strerror(errno)) << " "
#define log_error std::cerr << "[ERROR] (" << __FILE__ << ": " << __LINE__ << ":errno: " << clean_errno()

  constexpr uint8_t FLXWordLength = 32;

  template < typename A, typename B >
  bool comp(A a, B b)
  {
    return a.second < b.second;
  }
} //namespace mvtx_utils

#endif
