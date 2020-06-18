// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_TDIRECTORYHELPER_H
#define FUN4ALL_TDIRECTORYHELPER_H

#include <string>
#include <vector>

class TH1;
class TDirectory;
class TFile;

class TDirectoryHelper
{
 public:
  static void duplicateDir(TDirectory* dest, TDirectory* source);

  static TH1* getHisto(TDirectory* dir, const std::string& histoname,
                       const std::string& where);

  static TDirectory* mkdir(TDirectory* topDir,
                           const std::string& path,
                           std::vector<std::string>* titles = 0);

  static bool mkpath(TDirectory* dir, const std::string& path);

  static void copyToFile(TDirectory* src, TFile* dest);

  static bool pathIsInDir(const std::string& path, TDirectory* dir);

  static void splitPath(const std::string& path,
                        std::vector<std::string>& paths);
};

#endif
