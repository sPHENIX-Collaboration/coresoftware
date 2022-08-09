// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
// -- Purpose: Test the same number of events written are read by protobufIO
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
  for (int i = 0; i < 123456; ++i) {
    evt.set_event_number(i);
    writer->write_event(evt);
  }
  writer->close();

  std::cout << "Written 123456 empty events to a stringstream, sized: "
            << (binstream.str().size() / (1024*1024)) << " MB." << std::endl;

  auto reader = std::make_shared<HepMC3::Readerprotobuf>(binstream);
  HepMC3::GenEvent evt_in;
  for (int i = 0; i < 123456; ++i) {
    if (!reader->read_event(evt_in)) {
      std::cout << "[ERROR]: Failed to read event " << i << std::endl;
      return 1;
    }
    if (evt_in.event_number() != i) {
      std::cout << "[ERROR]: Read event " << i
                << ", but it had event number: " << evt_in.event_number()
                << std::endl;
      return 2;
    }
  }

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

  std::stringstream ss1("");
  std::stringstream ss2("");

  HepMC3::Print::listing(ss1, evt_in);
  HepMC3::Print::listing(ss2, evt);
  std::cout << ss1.str() << "\n" << ss2.str() << std::endl;
  if (ss1.str() != ss2.str()) {
    return 1;
  }

  return 0;
}
