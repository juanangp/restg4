

#include "TrackingAction.h"

#include "SimulationManager.h"

using namespace std;

TrackingAction::TrackingAction(SimulationManager* simulationManager)
    : G4UserTrackingAction(), fSimulationManager(simulationManager) {}

TrackingAction::~TrackingAction() {}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {
    fSimulationManager->GetOutputManager()->RecordTrack(track);
}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
    fSimulationManager->GetOutputManager()->UpdateTrack(track);
}
