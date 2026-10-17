//
// Created by lobis on 6/20/2022.
//

#include "SteppingVerbose.h"

#include <SimulationManager.h>



SteppingVerbose::SteppingVerbose(SimulationManager* simulationManager)
    : G4SteppingVerbose(), fSimulationManager(simulationManager) {}
SteppingVerbose::~SteppingVerbose() {}

void SteppingVerbose::TrackingStarted() {
    CopyState();
    fSimulationManager->GetOutputManager()->RecordStep(fStep);
}

void SteppingVerbose::StepInfo() {}
