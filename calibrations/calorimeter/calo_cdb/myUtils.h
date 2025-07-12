#ifndef CALOCDB_MYUTILS_H
#define CALOCDB_MYUTILS_H

// ROOT includes --
#include <TFitResultPtr.h>
#include <TH1.h>

// -- c++ includes --
#include <concepts>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

template <typename Func>
concept InvocableWithString = std::invocable<Func, const std::string&>;

class myUtils
{
 public:
  static void setEMCalDim(TH1* hist);

  static std::pair<int, int> getSectorIB(int iphi, int ieta);
  static std::pair<int, int> getSectorIB(int towerIndex);
  static std::vector<std::string> split(const std::string& s, char delimiter);
  static TFitResultPtr doGausFit(TH1* hist, Double_t start, Double_t end, const std::string& name = "fitFunc");

  /**
   * @brief Reads a CSV (or any line-delimited) file and applies a handler function to each line.
   *
   * @tparam Callable The type of the function/lambda to be called for each line.
   * Must be invocable with a 'const std::string&'.
   * @param filePath The path to the input file.
   * @param lineHandler A function, lambda, or functor that takes a 'const std::string&' (the line)
   * and processes it.
   * @param skipHeader If true, the first line of the file will be read and discarded.
   * @return true if the file was successfully opened and read, false otherwise.
   */
  template <InvocableWithString Callable>  // Using the more general concept for wider applicability
  static Bool_t readCSV(const std::filesystem::path& filePath, Callable lineHandler, Bool_t skipHeader = true)
  {
    std::ifstream file(filePath);

    if (!file.is_open())
    {
      std::cout << "Error: [" << filePath.string() << "] Could not open file." << std::endl;
      return false;
    }

    std::string line;

    if (skipHeader && std::getline(file, line))
    {
      // First line read and discarded (header)
    }

    while (std::getline(file, line))
    {
      // Optional: Handle potential Windows CRLF (\r\n) issues if the file might
      // have them and you're on a system that only expects \n.
      // std::getline usually handles this, but if \r remains:
      // if (!line.empty() && line.back() == '\r') {
      //     line.pop_back();
      // }
      lineHandler(line);  // Call the user-provided function for each line
    }

    // Check for errors during read operations (other than EOF)
    if (file.bad())
    {
      std::cout << "Error: [" << filePath.string() << "] I/O error while reading file." << std::endl;
      return false;
    }
    // file.eof() will be true if EOF was reached.
    // file.fail() might be true if getline failed not due to eof (e.g., badbit also set).
    // If badbit is not set and eof() is true, it's a successful read to the end.

    return true;  // Successfully processed or reached EOF
  }
};

#endif
