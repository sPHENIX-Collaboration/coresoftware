// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_COMPRESSEDIO_H
#define HEPMC3_COMPRESSEDIO_H
#if HEPMC3_USE_COMPRESSION
#if HEPMC3_Z_SUPPORT
#define BXZSTR_Z_SUPPORT 1
#endif
#if HEPMC3_LZMA_SUPPORT
#define BXZSTR_LZMA_SUPPORT 1
#endif
#if HEPMC3_BZ2_SUPPORT
#define BXZSTR_BZ2_SUPPORT 1
#endif
#include "HepMC3/bxzstr/bxzstr.hpp"
namespace HepMC3
{
using ofstream = bxz::ofstream;
using ostream = bxz::ostream;
using ifstream = bxz::ifstream;
using istream = bxz::istream;

using Compression = bxz::Compression;
inline Compression detect_compression_type(char* in_buff_start, char* in_buff_end) {
    return bxz::detect_type(in_buff_start,in_buff_end);
}
const  std::vector<Compression> supported_compression_types = {
#if HEPMC3_Z_SUPPORT
    Compression::z,
#endif
#if HEPMC3_LZMA_SUPPORT
    Compression::lzma,
#endif
#if HEPMC3_BZ2_SUPPORT
    Compression::bz2,
#endif
};
std::vector<Compression> known_compression_types = {
    Compression::z,
    Compression::lzma,
    Compression::bz2
};
}
namespace std
{
string to_string(HepMC3::Compression & c) {
    switch (c) {
    case HepMC3::Compression::z:
        return string("z");
    case HepMC3::Compression::lzma:
        return string("lzma");
    case HepMC3::Compression::bz2:
        return string("bz2");
    default:
        break;
    }
    return string("plaintext");
}
}

#endif
#endif
