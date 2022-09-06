#include <iostream>
#include "AliAnalysisPairExtractor.h"
#include "AliAnalysisPairExtractor.cxx"

using std::cout;
using std::endl;
using std::flush;

void PairTreeMakerLaptop() {

  cout << "Beginn!" << endl;
  
  AliAnalysisPairExtractor ana;
  ana.setPDG(443, -11, 11);
  ana.SetUp("tree_candidates.root");
  cout << "Object Created!" << endl;
  
  ana.extractDirectory("/alice/data/runs", "JpsiCandidates_data.root", "DstTree", "MB", "Minimum Bias Events Flag 14.", kFALSE);
  cout << "data extraction finished!" << endl;

  ana.extractDirectory("/alice/sim/runs", "JpsiCandidates_MC.root", "DstTree", "Injected", "injected Monte Carlo.", kTRUE);
  cout << "mc extraction finished!" << endl;
  
  ana.Write();
  cout << "File saved!" << endl;
}
