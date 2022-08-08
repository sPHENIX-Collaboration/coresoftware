// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test that we can correctly deduce the reader type when passed a
// binary protobuf file
//

// These are the only headers in ReaderPlugin, so including these firstmakes
// sure the hack below doesn't leak out of the ReaderPlugin header
#include "HepMC3/Reader.h"
#include "HepMC3/ReaderFactory.h"
int main() {
  std::shared_ptr<HepMC3::Reader> rdr = HepMC3::deduce_reader("inputIO20.proto");
  if (!rdr) return 1;
  rdr->close();
  return 0;
}
