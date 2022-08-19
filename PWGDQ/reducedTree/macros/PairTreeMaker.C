#include <iostream>
#include <array>
#include <string>
#include <dirent.h>
#include <AliAnalysisPairExtractor.h>
#include "AliAnalysisPairExtractor.cxx"

using std::cout;
using std::endl;
using std::flush;

void PairTreeMaker() {
	
  cout << "Beginn!" << endl;
  
  AliAnalysisPairExtractor ana;
  ana.setPDG(443, -11, 11);
  ana.SetUp("tree_candidates.root");
    
  cout << "Object Created!" << endl;
  
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir("/alice/data/runs")) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir(dir)) != NULL) {
      string name(ent->d_name);
      if (name.find(".") != std::string::npos) continue;
      ana.extractDataFile("/alice/data/runs/" + name + "/JpsiCandidates_data.root");
      cout << "Data Run: " << name << endl;
    }
    closedir (dir);
  }
  if ((dir = opendir("/alice/sim/runs")) != NULL) {
    /* print all the files and directories within directory */
    while ((ent = readdir(dir)) != NULL) {
      string name(ent->d_name);
      if (name.find(".") != std::string::npos) continue;
      ana.extractMCFile("/alice/sim/runs/" + name + "/JpsiCandidates_MC.root");
      cout << "MC Run: " << name << endl;
    }
    closedir (dir);
  }

  cout << "data extraction finished!" << endl;
  
  ana.Write();
  
  cout << "File saved!" << endl;
}
