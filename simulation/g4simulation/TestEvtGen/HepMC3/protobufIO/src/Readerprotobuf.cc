// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
/**
 *  @file Readerprotobuf.cc
 *  @brief Implementation of \b class Readerprotobuf
 *
 */
#include "HepMC3/Readerprotobuf.h"

#include "HepMC3/Print.h"
#include "HepMC3/Version.h"

#include "HepMC3/Data/GenRunInfoData.h"

// protobuf header files
#include "HepMC3.pb.h"

namespace HepMC3 {

std::string const ProtobufMagicHeader = "hmpb";
size_t const ProtobufMagicHeaderBytes = 4;

HEPMC3_DECLARE_READER_FILE(Readerprotobuf);
HEPMC3_DECLARE_READER_STREAM(Readerprotobuf);

static size_t const MDBytesLength = 10;

Readerprotobuf::Readerprotobuf(const std::string &filename)
    : m_in_file(nullptr), m_bytes_read(0),
      m_msg_type(HepMC3_pb::MessageDigest::unknown) {

  m_md_buffer.resize(MDBytesLength);

  m_in_file = std::unique_ptr<std::ifstream>(
      new std::ifstream(filename, ios::in | ios::binary));

  if (!m_in_file->is_open()) {
    HEPMC3_ERROR("Readerprotobuf: Problem opening file: " << filename)
    return;
  }

  m_in_stream = m_in_file.get();

  read_file_start();
}

Readerprotobuf::Readerprotobuf(std::istream &stream)
    : m_in_file(nullptr), m_bytes_read(0),
      m_msg_type(HepMC3_pb::MessageDigest::unknown) {

  if (!stream.good()) {
    HEPMC3_ERROR(
        "Cannot initialize Readerprotobuf on istream which is not good().");
    return;
  }

  m_md_buffer.resize(MDBytesLength);

  m_in_stream = &stream;
  read_file_start();
}

Readerprotobuf::Readerprotobuf(std::shared_ptr<std::istream> stream)
    : Readerprotobuf(*stream) {}

bool Readerprotobuf::read_file_start() {

  // Read the first 16 bytes, it should read "HepMC3::Protobuf"
  std::string MagicIntro;
  MagicIntro.resize(ProtobufMagicHeaderBytes);
  m_in_stream->read(&MagicIntro[0], ProtobufMagicHeaderBytes);

  if (MagicIntro != ProtobufMagicHeader) {
    HEPMC3_ERROR("Failed to find expected Magic first "
                 << ProtobufMagicHeaderBytes
                 << " bytes, is this really "
                    "a HepMC3::Protobuf file?");
    return false;
  }

  if (!read_Header()) {
    HEPMC3_ERROR("Readerprotobuf: Problem parsing start of file, expected to "
                 "find Header, but instead found message type: "
                 << m_msg_type);
    return false;
  }

  if (!read_GenRunInfo()) {
    HEPMC3_ERROR("Readerprotobuf: Problem parsing start of file, expected to "
                 "find RunInfo, but instead found message type: "
                 << m_msg_type);
    return false;
  }

  return true;
}

bool Readerprotobuf::buffer_message() {
  if (failed()) {
    return false;
  }

  if (m_msg_buffer.size()) { // if we already have a message that hasn't been
                             // parsed, don't buffer the next one
    return true;
  }

  m_msg_type = HepMC3_pb::MessageDigest::unknown;

  m_in_stream->read(&m_md_buffer[0], MDBytesLength);

  if (failed()) {
    return false;
  }

  m_bytes_read += MDBytesLength;

  HepMC3_pb::MessageDigest md;
  if (!md.ParseFromString(m_md_buffer)) {
    return false;
  }

  m_msg_type = md.message_type();

  m_msg_buffer.resize(md.bytes());
  m_in_stream->read(&m_msg_buffer[0], md.bytes());

  if (failed()) {
    return false;
  }

  m_bytes_read += md.bytes();

  if (m_msg_type ==
      HepMC3_pb::MessageDigest::Footer) { // close the stream if we have read to
                                          // the end of the file
    close();
  }

  return true;
}

bool Readerprotobuf::read_Header() {
  if (!buffer_message()) {
    return false;
  }

  if (m_msg_type != HepMC3_pb::MessageDigest::Header) {
    return false;
  }

  HepMC3_pb::Header Header_pb;
  if (!Header_pb.ParseFromString(m_msg_buffer)) {
    // if we fail to read a message then close the stream to indicate failed
    // state
    close();
    return false;
  }
  m_msg_buffer.clear();

  m_file_header.m_version_str = Header_pb.version_str();
  m_file_header.m_version_maj = Header_pb.version_maj();
  m_file_header.m_version_min = Header_pb.version_min();
  m_file_header.m_version_patch = Header_pb.version_patch();
  m_file_header.m_protobuf_version_maj = Header_pb.protobuf_version_maj();
  m_file_header.m_protobuf_version_min = Header_pb.protobuf_version_min();
  m_file_header.m_protobuf_version_patch = Header_pb.protobuf_version_patch();

  return true;
}

bool Readerprotobuf::read_GenRunInfo() {
  if (!buffer_message()) {
    return false;
  }

  if (m_msg_type != HepMC3_pb::MessageDigest::RunInfo) {
    return false;
  }

  set_run_info(std::make_shared<HepMC3::GenRunInfo>());

  HepMC3_pb::GenRunInfoData GenRunInfo_pb;
  if (!GenRunInfo_pb.ParseFromString(m_msg_buffer)) {
    // if we fail to read a message then close the stream to indicate failed
    // state
    close();
    return false;
  }
  m_msg_buffer.clear();

  HepMC3::GenRunInfoData gridata;

  int vector_size = 0;

  vector_size = GenRunInfo_pb.weight_names_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.weight_names.push_back(GenRunInfo_pb.weight_names(it));
  }

  vector_size = GenRunInfo_pb.tool_name_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.tool_name.push_back(GenRunInfo_pb.tool_name(it));
  }

  vector_size = GenRunInfo_pb.tool_version_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.tool_version.push_back(GenRunInfo_pb.tool_version(it));
  }

  vector_size = GenRunInfo_pb.tool_description_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.tool_description.push_back(GenRunInfo_pb.tool_description(it));
  }

  vector_size = GenRunInfo_pb.attribute_name_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.attribute_name.push_back(GenRunInfo_pb.attribute_name(it));
  }

  vector_size = GenRunInfo_pb.attribute_string_size();
  for (int it = 0; it < vector_size; ++it) {
    gridata.attribute_string.push_back(GenRunInfo_pb.attribute_string(it));
  }

  run_info()->read_data(gridata);
  return true;
}

bool Readerprotobuf::read_GenEvent(bool skip) {
  if (!buffer_message()) {
    return false;
  }

  if (m_msg_type != HepMC3_pb::MessageDigest::Event) {
    return false;
  }

  if (skip) { // Don't parse to HepMC3 if skipping
    m_msg_buffer.clear();
    return true;
  }

  if (!m_msg_buffer.size()) { // empty event
    m_evdata = HepMC3::GenEventData();
    return true;
  }

  HepMC3_pb::GenEventData ged_pb;
  if (!ged_pb.ParseFromString(m_msg_buffer)) {
    // if we fail to read a message then close the stream to indicate failed
    // state
    close();
    return false;
  }

  m_evdata.event_number = ged_pb.event_number();

  switch (ged_pb.momentum_unit()) {
  case HepMC3_pb::GenEventData::MEV: {
    m_evdata.momentum_unit = HepMC3::Units::MEV;
    break;
  }
  case HepMC3_pb::GenEventData::GEV: {
    m_evdata.momentum_unit = HepMC3::Units::GEV;
    break;
  }
  default: {
    HEPMC3_ERROR("Unknown momentum unit: " << ged_pb.momentum_unit());
    return false;
  }
  }

  switch (ged_pb.length_unit()) {
  case HepMC3_pb::GenEventData::MM: {
    m_evdata.length_unit = HepMC3::Units::MM;
    break;
  }
  case HepMC3_pb::GenEventData::CM: {
    m_evdata.length_unit = HepMC3::Units::CM;
    break;
  }
  default: {
    HEPMC3_ERROR("Unknown length unit: " << ged_pb.length_unit());
    return false;
  }
  }

  int vector_size = 0;

  m_evdata.particles.clear();
  vector_size = ged_pb.particles_size();
  for (int it = 0; it < vector_size; ++it) {
    auto particle_pb = ged_pb.particles(it);

    HepMC3::GenParticleData pdata;

    pdata.pid = particle_pb.pid();
    pdata.status = particle_pb.status();
    pdata.is_mass_set = particle_pb.is_mass_set();
    pdata.mass = particle_pb.mass();

    pdata.momentum = HepMC3::FourVector{
        particle_pb.momentum().m_v1(), particle_pb.momentum().m_v2(),
        particle_pb.momentum().m_v3(), particle_pb.momentum().m_v4()};

    m_evdata.particles.push_back(pdata);
  }

  m_evdata.vertices.clear();
  vector_size = ged_pb.vertices_size();
  for (int it = 0; it < vector_size; ++it) {
    auto vertex_pb = ged_pb.vertices(it);

    HepMC3::GenVertexData vdata;

    vdata.status = vertex_pb.status();

    vdata.position = HepMC3::FourVector{
        vertex_pb.position().m_v1(), vertex_pb.position().m_v2(),
        vertex_pb.position().m_v3(), vertex_pb.position().m_v4()};

    m_evdata.vertices.push_back(vdata);
  }

  m_evdata.weights.clear();
  vector_size = ged_pb.weights_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.weights.push_back(ged_pb.weights(it));
  }

  m_evdata.links1.clear();
  vector_size = ged_pb.links1_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.links1.push_back(ged_pb.links1(it));
  }

  m_evdata.links2.clear();
  vector_size = ged_pb.links2_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.links2.push_back(ged_pb.links2(it));
  }

  m_evdata.event_pos =
      HepMC3::FourVector{ged_pb.event_pos().m_v1(), ged_pb.event_pos().m_v2(),
                         ged_pb.event_pos().m_v3(), ged_pb.event_pos().m_v4()};

  m_evdata.attribute_id.clear();
  vector_size = ged_pb.attribute_id_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.attribute_id.push_back(ged_pb.attribute_id(it));
  }

  m_evdata.attribute_name.clear();
  vector_size = ged_pb.attribute_name_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.attribute_name.push_back(ged_pb.attribute_name(it));
  }

  m_evdata.attribute_string.clear();
  vector_size = ged_pb.attribute_string_size();
  for (int it = 0; it < vector_size; ++it) {
    m_evdata.attribute_string.push_back(ged_pb.attribute_string(it));
  }

  m_msg_buffer.clear();
  return true;
}

bool Readerprotobuf::skip(const int n) {

  for (int nn = n; nn > 0; --nn) {
    if (!read_GenEvent(true)) {
      return false;
    }
  }
  return !failed();
}

bool Readerprotobuf::read_event(GenEvent &evt) {

  if (!read_GenEvent(false)) {
    return false;
  }

  evt.read_data(m_evdata);
  evt.set_run_info(run_info());

  return true;
}

void Readerprotobuf::close() {
  if (m_in_file) {
    m_in_file->close();
    m_in_file.reset();
  }
  m_in_stream = nullptr;
  m_msg_buffer.clear();
}

bool Readerprotobuf::failed() {
  if (m_in_file) {
    return !m_in_file->is_open() || !m_in_file->good();
  }
  return !m_in_stream || !m_in_stream->good();
}

} // namespace HepMC3
