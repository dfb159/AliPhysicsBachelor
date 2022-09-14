#include <iostream>
#include <AliAnalysisPairExtractor.h>

using std::cout;
using std::endl;
using std::flush;

void PairTreeMakerLaptop() {
  AliAnalysisPairExtractor ana("tree_candidates.root");
  cout << "Object Created!" << endl;

  ana.setOutputTree("data_MB", "Minimum Bias GSI data with Event Flag 14");
  ana.extractDirectory("/alice/data/runs", "JpsiCandidates_data.root", "DstTree");
  cout << "data extraction finished!" << endl;
  
  ana.addPDGCut(443, -11, 11, kTRUE);
  ana.setOutputTree("mc_Injected_true", "injected Monte Carlo with confirmed origin");
  ana.extractDirectory("/alice/sim/runs", "JpsiCandidates_MC.root", "DstTree");
  ana.clearFilters();
  cout << "mc true extraction finished!" << endl;
  
  ana.addPDGCut(443, -11, 11, kFALSE);
  ana.setOutputTree("mc_Injected_false", "injected Monte Carlo with rejected origin");
  ana.extractDirectory("/alice/sim/runs", "JpsiCandidates_MC.root", "DstTree");
  ana.clearFilters();
  cout << "mc true extraction finished!" << endl;
}
