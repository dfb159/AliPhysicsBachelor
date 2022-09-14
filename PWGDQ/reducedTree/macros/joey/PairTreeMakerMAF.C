#include <iostream>
#include <AliAnalysisPairExtractor.h>

using std::cout;
using std::endl;
using std::flush;

void PairTreeMakerMAF() {
  AliAnalysisPairExtractor ana("/gluster1/j_sigr01/analysis/tree_candidates.root");
  cout << "Analysis created!" << endl;
  
  ana.setOutputTree("data_MB", "Minimum Bias GSI data with Event Flag 14");
  ana.extractDirectory("/gluster1/j_sigr01/data_gsi/runs", "JpsiCandidates_data.root", "DstTree");
  cout << "data mb extraction finished!" << endl;

  ana.setOutputTree("data_Enriched", "All enriched GSI data");
  ana.extractDirectory("/gluster1/j_sigr01/data_gsi_all/runs", "JpsiCandidates_data.root", "DstTree");
  cout << "data enriched extraction finished!" << endl;

  ana.addPDGCut(443, -11, 11, kTRUE);
  ana.setOutputTree("mc_Injected_true", "injected Monte Carlo with confirmed origin");
  ana.extractDirectory("/gluster1/j_sigr01/MC_injected/runs", "JpsiCandidates_MC.root", "DstTree");
  ana.clearFilters();
  cout << "mc true extraction finished!" << endl;
  
  ana.addPDGCut(443, -11, 11, kFALSE);
  ana.setOutputTree("mc_Injected_false", "injected Monte Carlo with rejected origin");
  ana.extractDirectory("/gluster1/j_sigr01/MC_injected/runs", "JpsiCandidates_MC.root", "DstTree");
  ana.clearFilters();
  cout << "mc false extraction finished!" << endl;
}
