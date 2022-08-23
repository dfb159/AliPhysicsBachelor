#include <iostream>
#include <AliAnalysisPairExtractor.h>
#include "AliAnalysisPairExtractor.cxx"

using std::cout;
using std::endl;
using std::flush;

void PairTreeMaker(TString outfile="tree_candidates.root", TString dataDir="/alice/data/runs", TString mcDir="/alice/sim/runs") {
	
  cout << "Beginn!" << endl;
  
  AliAnalysisPairExtractor ana;
  ana.setPDG(443, -11, 11);
  ana.SetUp(outfile);
    
  cout << "Object Created!" << endl;
  
  ana.extractDataDirectory(dataDir);
  cout << "data extraction finished!" << endl;

  ana.extractMCDirectory(mcDir);
  cout << "mc extraction finished!" << endl;
  
  ana.Write();
  
  cout << "File saved!" << endl;
}
