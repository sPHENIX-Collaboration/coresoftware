// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2022 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERPROTOBUF_H
#define HEPMC3_READERPROTOBUF_H
/**
 *  @file  Readerprotobuf.h
 *  @brief Definition of \b class Readerprotobuf
 *
 *  @class HepMC3::Readerprotobuf
 *  @brief GenEvent I/O parsing and serialization for protobuf-based binary
 * files
 *
 * If HepMC was compiled with path to protobuf available, this class can be
 * used for protobuf file I/O in the same manner as with HepMC::ReaderAscii
 * class.
 *
 *  @ingroup IO
 *
 */

#include "HepMC3/Reader.h"

#include "HepMC3/GenEvent.h"

#include "HepMC3/Data/GenEventData.h"

#include <array>
#include <fstream>
#include <string>
#include <vector>

namespace HepMC3 {

class Readerprotobuf : public Reader {
public:
  /**
   * @class HepMC3::Readerprotobuf::FileHeader
   * @brief A copy of the information contained in the protobuf file header
   */
  struct FileHeader {
    std::string m_version_str;
    unsigned int m_version_maj;
    unsigned int m_version_min;
    unsigned int m_version_patch;

    unsigned int m_protobuf_version_maj;
    unsigned int m_protobuf_version_min;
    unsigned int m_protobuf_version_patch;
  };

  //
  // Constructors
  //
public:
  /** @brief filename constructor
   *
   * @details Attempts to open the passed filename and read protobuf HepMC3
   * events from it
   */
  Readerprotobuf(const std::string &filename);

  /** @brief istream constructor
   *
   * @details Attempts to read a binary HepMC3 protobuf event stream from the
   * passed istream object
   */
  Readerprotobuf(std::istream &stream);

  /** @brief istream constructor
   *
   * @details Attempts to read a binary HepMC3 protobuf event stream from the
   * passed istream object
   */
  Readerprotobuf(std::shared_ptr<std::istream> stream);

  //
  // Functions
  //
public:
  /** @brief skips the next N events
   *
   *  @param[in] the number of events to skip
   *  @return Whether the reader can still be read from after skipping
   */
  bool skip(const int) override;

  /** @brief Read event from file
   *
   *  @param[out] evt Contains parsed event
   *  @return Whether the reader can still be read from after reading
   */
  bool read_event(GenEvent &evt) override;

  /** @brief Close file stream */
  void close() override;

  /** @brief Get the header information read from the protobuf file */
  FileHeader const &file_header() { return m_file_header; }

  /** @brief Get stream error state */
  bool failed() override;
  //
  // Fields
  //
private:
  /** @brief Read the next protobuf message into the message buffer
   *
   * @details Fills m_msg_buffer with the next message and sets m_msg_type to
   * signify the message type. Returns true if there is a message in the buffer
   * ready to parse.
   */
  bool buffer_message();

  /** @brief Parse the next protobuf message as a GenRunInfo message
   *
   * @return Whether the reader can still be read from after reading
   */
  bool read_GenRunInfo();

  /** @brief Parse the next protobuf message as a GenEvent message
   *
   * @param[in] skip Whether to bother actually parsing this message to a
   * GenEvent
   * @return Whether the reader can still be read from after reading
   */
  bool read_GenEvent(bool skip = false);

  /** @brief Parse the next protobuf message as a Header message
   *
   * @return Whether the reader can still be read from after reading
   */
  bool read_Header();

  /** @brief Parse the front matter of the protobuf message stream before the
   * events
   */
  bool read_file_start();

  /** @brief The total number of event bytes read, including message frames
   */
  size_t m_bytes_read;

  /** @brief The file stream of the file being read
   * 
   * @detail This is non-null and owned by this class if constructed with the 
   * string constructor, otherwise it will be null
   */
  std::unique_ptr<std::ifstream> m_in_file;
  /** @brief The stream object that is read from
   * 
   * @detail If constructed with the string constructor, this just points to 
   * m_in_file.get())
   */
  std::istream *m_in_stream;

  /** @brief The buffer used to hold the current message binary 
   * (header/genruninfo/genevent/footer)
   */
  std::string m_msg_buffer;
  /** @brief The buffer used to hold the current message digest binary 
   * (message frame)
   */
  std::string m_md_buffer;
  /** @brief The type of current message
   * 
   * @details Defined in HepMC3_pb::MessageDigest::MessageType in the proto file
   */
  int m_msg_type;

  /** @brief The event data parsed from the message
   */
  HepMC3::GenEventData m_evdata;

  /** @brief A copy of the library version info stored in the proto file header
   * 
   * @detail This is a copy so as to avoid passing on protobuf header 
   * dependencies to files that include this header
   */
  FileHeader m_file_header;
};

} // namespace HepMC3

#endif
