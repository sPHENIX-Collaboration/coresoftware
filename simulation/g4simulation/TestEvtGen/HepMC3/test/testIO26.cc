// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test that an event read from an ASCII input file can be
// serialized to and from protobuf and still be matched to the same event
//

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/Readerprotobuf.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/Writerprotobuf.h"
#include "HepMC3TestUtils.h"

using namespace HepMC3;

int main() {
  ReaderAsciiHepMC2 inputA("inputIO26.hepmc");
  if (inputA.failed()) {
    return 1;
  }

  Writerprotobuf outputA("frominputIO26.proto");
  if (outputA.failed()) {
    return 2;
  }

  while (!inputA.failed()) {
    GenEvent evt_in_ASCII(Units::GEV, Units::MM);
    inputA.read_event(evt_in_ASCII);
    if (inputA.failed()) {
      printf("End of file reached. Exit.\n");
      break;
    }
    outputA.write_event(evt_in_ASCII);
  }
  inputA.close();
  outputA.close();

  Readerprotobuf inputB("frominputIO26.proto");
  if (inputB.failed()) {
    return 3;
  }

  WriterAsciiHepMC2 outputB("fromfrominputIO26.hepmc");
  if (outputB.failed()) {
    return 4;
  }

  while (!inputB.failed()) {
    GenEvent evt_in_proto(Units::GEV, Units::MM);
    inputB.read_event(evt_in_proto);
    if (inputB.failed()) {
      printf("End of file reached. Exit.\n");
      break;
    }
    outputB.write_event(evt_in_proto);
  }

  inputB.close();
  outputB.close();

  return COMPARE_ASCII_FILES("fromfrominputIO26.hepmc", "inputIO26.hepmc");
}
