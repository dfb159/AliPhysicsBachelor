#include <iostream>
#include "AliAnalysisPairExtractor.h"
#include "AliAnalysisPairExtractor.cxx"

using std::cout;
using std::endl;
using std::flush;

void PairTreeMakerMAF() {

  cout << "Beginn!" << endl;
  
  AliAnalysisPairExtractor ana;
  ana.setPDG(443, -11, 11);
  ana.SetUp("/gluster1/j_sigr01/analysis/tree_candidates.root");
  cout << "Object Created!" << endl;
  
  ana.extractDirectory("/gluster1/j_sigr01/data_gsi/runs", "JpsiCandidates_data.root", "DstTree", "MB", "Minimum Bias GSI data with Event Flag 14.", kFALSE);
  cout << "data mb extraction finished!" << endl;

  ana.extractDirectory("/gluster1/j_sigr01/data_gsi_all/runs", "JpsiCandidates_data.root", "DstTree", "Enriched", "All enriched GSI data", kFALSE);
  cout << "data enriched extraction finished!" << endl;

  ana.extractDirectory("/gluster1/j_sigr01/MC_injected/runs", "JpsiCandidates_MC.root", "DstTree", "Injected", "injected Monte Carlo.", kTRUE);
  cout << "mc extraction finished!" << endl;
  
  ana.Write();
  cout << "File saved!" << endl;
}
