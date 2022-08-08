// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test graceful failure for impossible-to-open file for writer
//

#include "HepMC3/Writerprotobuf.h"

int main() {
  HepMC3::Writerprotobuf wrtr(
      "/path/to/invalid/place/erihreiuhgeihgir/test.proto");
  if (!wrtr.failed()) return 1;
  return 0;
}
