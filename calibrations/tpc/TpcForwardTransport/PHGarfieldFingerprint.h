#ifndef PHGARFIELDFINGERPRINT_H
#define PHGARFIELDFINGERPRINT_H

#include <bit>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>

// Internal byte-order contract shared by the existing map fingerprints.
namespace PHGarfieldFingerprint
{
  class Fnv1a64
  {
   public:
    void Byte(unsigned char value)
    {
      m_value ^= value;
      m_value *= 1099511628211ULL;
    }
    void Integer(std::uint64_t value)
    {
      for (unsigned int i = 0; i < 8; ++i)
      {
        Byte(static_cast<unsigned char>(value & 255U));
        value >>= 8;
      }
    }
    void Real(double value) { Integer(std::bit_cast<std::uint64_t>(value)); }
    void Bytes(std::string_view value)
    {
      for (unsigned char character : value)
      {
        Byte(character);
      }
    }
    std::string Hex() const
    {
      std::ostringstream result;
      result << std::hex << std::setfill('0') << std::setw(16) << m_value;
      return result.str();
    }

   private:
    std::uint64_t m_value{14695981039346656037ULL};
  };
}  // namespace PHGarfieldFingerprint

#endif
