
#include "SteppingAction.h"

#include "SimulationManager.h"

using namespace std;

SteppingAction::SteppingAction(SimulationManager* simulationManager)
    : fSimulationManager(simulationManager) {}

SteppingAction::~SteppingAction() {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    fSimulationManager->GetOutputManager()->RecordStep(step);
}
