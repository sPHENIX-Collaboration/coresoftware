// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test the fail state of the reader and the return value of
// read_event near the end of the event stream
//

#include "HepMC3/Readerprotobuf.h"
#include "HepMC3/Writerprotobuf.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"

int main() {

  HepMC3::GenEvent evt;

  std::stringstream binstream;

  auto writer = std::make_shared<HepMC3::Writerprotobuf>(binstream);
  writer->write_event(evt);
  writer->write_event(evt);
  writer->write_event(evt);
  writer->close();

  auto reader = std::make_shared<HepMC3::Readerprotobuf>(binstream);
  HepMC3::GenEvent evt_in;
  reader->read_event(evt_in);
  reader->read_event(evt_in);
  reader->read_event(evt_in);

  // we have now read all events from file, ensure that the stream is still
  // considered good
  if (reader->failed()) {
    return 3;
  }
  // the next event read should fail
  if (reader->read_event(evt_in)) {
    return 4;
  }
  // and should set the stream to bad
  if (!reader->failed()) {
    return 5;
  }

  reader->close();
  return 0;
}
