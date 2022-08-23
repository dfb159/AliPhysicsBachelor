#include <iostream>
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
  
  ana.extractDataDirectory("/alice/data/runs");
  ana.extractMCDirectory("/alice/sim/runs");

  cout << "data extraction finished!" << endl;
  
  ana.Write();
  
  cout << "File saved!" << endl;
}
