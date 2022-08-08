// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_VERSION_H
#define HEPMC3_VERSION_H

#include <string>

/// HepMC version string
#define HEPMC3_VERSION "3.02.05"

/// @brief HepMC version as an integer, HepMC X.Y.Z = 1000000*X + 1000*Y + Z
///
/// Use like "#if HEPMC3_VERSION_CODE < 3001004" for < 3.01.04
#define HEPMC3_VERSION_CODE 3002005
namespace HepMC3 {
/// Get the HepMC library version string
inline std::string version() {
    return HEPMC3_VERSION;
}
}

#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
#define HEPMC3_DECLARE_READER_FILE(classname)    extern "C" {  __declspec(dllexport) classname *  __stdcall new ## classname ## file (const std::string &filename ) { return new classname (filename);  } }
#define HEPMC3_DECLARE_READER_STREAM(classname)  extern "C" {  __declspec(dllexport) classname *  __stdcall new ## classname ## stream (std::istream & stream) { return new classname (stream);  } }
#define HEPMC3_DECLARE_WRITER_FILE(classname)    extern "C" {  __declspec(dllexport) classname *  __stdcall new ## classname ## file (const std::string &filename, std::shared_ptr<GenRunInfo> run ) { return new classname (filename,run);  } }
#define HEPMC3_DECLARE_WRITER_STREAM(classname)  extern "C" {  __declspec(dllexport) classname * __stdcall new ## classname ## stream (std::ostream & stream, std::shared_ptr<GenRunInfo> run) { return new classname (stream,run);  } }
#endif
#if defined(__linux__) || defined(__darwin__)|| defined(__APPLE__) || defined(__FreeBSD__) || defined(__sun)
#define HEPMC3_DECLARE_READER_FILE(classname)    extern "C" { classname * new ## classname ## file (const std::string &filename ) { return new classname (filename);  } }
#define HEPMC3_DECLARE_READER_STREAM(classname)  extern "C" { classname * new ## classname ## stream (std::istream & stream) { return new classname (stream);  } }
#define HEPMC3_DECLARE_WRITER_FILE(classname)    extern "C" { classname * new ## classname ## file (const std::string &filename, std::shared_ptr<GenRunInfo> run ) { return new classname (filename,run);  } }
#define HEPMC3_DECLARE_WRITER_STREAM(classname)  extern "C" { classname * new ## classname ## stream (std::ostream & stream, std::shared_ptr<GenRunInfo> run) { return new classname (stream,run);  } }
#endif
#endif
