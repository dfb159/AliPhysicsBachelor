///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskFlatenicity macro                           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskFlatenicity *
AddTaskFlatenicity(const Char_t *taskname = "Flat", TString detForFlat = "V0",
                   Bool_t woTrivialscaling = kFALSE, Bool_t useMC = kTRUE,
                   Double_t minpT = 0.15)

{
  // detForFlat: "V0", "TPC", "V0_TPC"
  // get the manager via the static access member. since it's static, you don't
  // need an instance of the class to call the function

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and
  // feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // now you create an instance of your task
  AliAnalysisTaskFlatenicity *taskFlat =
      new AliAnalysisTaskFlatenicity("taskFlat");
  if (!taskFlat)
    return 0x0;
  taskFlat->SetUseMC(useMC);
  taskFlat->SetDetectorForFlatenicity(detForFlat);
  taskFlat->SetPtMin(minpT);
  taskFlat->SetRemoveTrivialScaling(woTrivialscaling);
  mgr->AddTask(taskFlat);

  const char *complement;
  if (woTrivialscaling) {
    complement = "wotrivialscal";
  } else {
    complement = "wtrivialscal";
  }
  const char *complement2;
  if (detForFlat == "V0") {
    complement2 = "V0";
  }
  if (detForFlat == "TPC") {
    complement2 = "TPC";
  }
  if (detForFlat == "V0_TPC") {
    complement2 = "V0_TPC";
  }

  // complement += detForFlat;

  mgr->ConnectInput(taskFlat, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
      taskFlat, 1,
      mgr->CreateContainer(
          Form("cList%s_%s_%s", taskname, complement, complement2),
          TList::Class(), AliAnalysisManager::kOutputContainer,
          Form("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname)));

  return taskFlat;
}
