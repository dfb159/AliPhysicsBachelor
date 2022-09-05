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
  
  ana.extractDirectory(dataDir, "JCandidates_data.root", "DstTree", "MB", "Minimum Bias Events Flag 14.", kFALSE);
  cout << "data extraction finished!" << endl;

  ana.extractDirectory(mcDir, "JCandidates_MC.root", "DstTree", "injected", "injected Monte Carlo.", kTRUE);
  cout << "mc extraction finished!" << endl;
  
  ana.Write();
  
  cout << "File saved!" << endl;
}
