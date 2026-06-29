#include <Garfield/MediumMagboltz.hh>
#include <phool/phool.h>

#include <filesystem>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <string>
#include <utility>
#include <sstream>
#include <iomanip>
#include <exception>


namespace fs = std::filesystem;
std::string mergedName(const std::string& path, unsigned int eindex);

bool searchAndUnpackDirectory(
    const std::string& directoryPath,
    std::set<unsigned int>& Eindices,
    std::set<unsigned int>& Bindices,
    std::map<std::pair<unsigned int, unsigned int>, std::string>& FileList);

int main(int argc, char* argv[])
{
  try
    {
      if (argc != 3)
	{
	  std::cerr << "Usage:\n"
		    << argv[0]
		    << " path_to_gasfiles name_of_output_file\n";
	  return 1;
	}
      
      const std::string path = argv[1];
      const std::string output = path + "/" + argv[2];
      
      std::set<unsigned int> Eindices;
      std::set<unsigned int> Bindices;
      std::map<std::pair<unsigned int, unsigned int>, std::string> FileList;
      
      if (!searchAndUnpackDirectory(path, Eindices, Bindices, FileList))
	{
	  std::cerr << PHWHERE << " Imperfect directory." << std::endl;
	  return 1;
	}
      
      Garfield::MediumMagboltz gas;
      
      for (const auto Eindex : Eindices)
	{
	  bool firstB = true;
	  
	  for (const auto Bindex : Bindices)
	    {
	      const auto it = FileList.find({Eindex, Bindex});
	      if (it == FileList.end())
		{
		  std::cerr << PHWHERE << " Missing file for E=" << Eindex
			    << " B=" << Bindex << std::endl;
		  return 1;
		}
	      
	      const std::string& nextfile = it->second;
	      
	      if (firstB)
		{
		  gas.LoadGasFile(nextfile);
		  firstB = false;
		}
	      else
		{
		  gas.MergeGasFile(nextfile, true);
		}
	    }
	  
	  const std::string mergedFile = mergedName(path, Eindex);
	  
	  std::cout << "Writing " << mergedFile << std::endl;
	  gas.WriteGasFile(mergedFile);
	  
	  std::vector<double> nE;
	  std::vector<double> nB;
	  std::vector<double> nA;
	  gas.GetFieldGrid(nE, nB, nA);
	  
	  std::cout << "Merged Gas File created: "<< mergedFile
		    << " with Grid Dimensions: " 
		    << nE.size() << " E-fields, " 
		    << nB.size() << " B-fields, " 
		    << nA.size() << " Angles." << std::endl;
	}
      
      bool firstE = true;
      
      for (const auto Eindex : Eindices)
	{
	  //if (Eindex > 10) {break;}
	  const std::string mergedFile = mergedName(path, Eindex);
	  
	  if (!fs::exists(mergedFile))
	    {
	      std::cerr << PHWHERE << " Missing merged file " << mergedFile << std::endl;
	      return 1;
	    }
	  
	  if (firstE)
	    {
	      gas.LoadGasFile(mergedFile);
	      firstE = false;
	    }
	  else
	    {
	      gas.MergeGasFile(mergedFile, true);
	    }
	}
      
      std::cout << "Writing final file " << output << std::endl;
      gas.WriteGasFile(output);
      
      std::vector<double> nE;
      std::vector<double> nB;
      std::vector<double> nA;
      gas.GetFieldGrid(nE, nB, nA);
      
      std::cout << "Final Gas File created: "<< output
		<< " with Grid Dimensions: " 
		<< nE.size() << " E-fields, " 
		<< nB.size() << " B-fields, " 
		<< nA.size() << " Angles." << std::endl;
      
      return 0;
    }

  catch (const std::exception& e)
    {
      std::cerr << PHWHERE << " Exception: " << e.what() << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << PHWHERE << " Unknown exception." << std::endl;
      return 1;
    }
}

bool searchAndUnpackDirectory(const std::string& directoryPath, std::set<unsigned int> &Eindices, std::set<unsigned int> &Bindices, std::map<std::pair<unsigned int, unsigned int>, std::string> &FileList)
{
    // Check if the directory exists and is valid
    if (!fs::exists(directoryPath) || !fs::is_directory(directoryPath))
      {
        std::cerr << "Error: Invalid directory path." << std::endl;
        return false;
      }

    std::regex filePattern(R"(^E([0-9]{3})_B([0-9]{3})\.gas$)");
    std::smatch matchResults;

    // Iterate through all items in the directory
    for (const auto& entry : fs::directory_iterator(directoryPath))
      {
        // Only process regular files
        if (entry.is_regular_file())
	  {
	    std::string filename = entry.path().filename().string();
	    
	    // Check if the filename matches our target pattern
	    if (std::regex_match(filename, matchResults, filePattern))
	      {
		// matchResults[1] contains the string after 'E'
		// matchResults[2] contains the string after 'B'
		// std::stoul automatically handles leading zeros
		unsigned int eValue = std::stoul(matchResults[1].str());
		unsigned int bValue = std::stoul(matchResults[2].str());
		Eindices.insert(eValue);
		Bindices.insert(bValue);
		FileList[{eValue, bValue}] = entry.path().string();
	      }
	  }
      }

    // Validate the results.
    if (Eindices.empty()) { return false; }
    if (Bindices.empty()) { return false; }

    unsigned int maxE = *Eindices.rbegin(); 
    unsigned int maxB = *Bindices.rbegin();
    for (unsigned int i=0; i<=maxE; i++)
      {
	for (unsigned int j=0; j<=maxB; j++)
	  {
	    if ( !FileList.contains({i,j}) ) { return false; }
	  }
      }

    std::cout << "   *** Gas File List Valid ***" << std::endl;
    std::cout << "Electric field indices 0 --> " << maxE << std::endl;
    std::cout << "Magnetic field indices 0 --> " << maxB << std::endl;
    
    return true;
}

std::string mergedName(const std::string& path, unsigned int eindex)
{
  std::ostringstream name;
  name << path << "/MERGED_E"
       << std::setw(3) << std::setfill('0') << eindex
       << ".gas";
  return name.str();
}
